#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/point-to-point-module.h"
#include "ns3/applications-module.h"
#include "ns3/ipv4-flow-classifier.h"
#include "ns3/flow-monitor-module.h"
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <limits>
#include <fstream>

using namespace ns3;

struct Edge {
    int src, dest;
};

// 定义个体
class DOA_Unit
{
public:
    std::vector<double> position; // 节点位置
    double value; // 适应度值
    double survival_rate;
    std::vector<std::vector<bool>> connections; // 节点连接矩阵
};

class DOA_Base
{
public:
    DOA_Base(int num_nodes, int size, int iter_max, const std::vector<double>& range_min_list, const std::vector<double>& range_max_list)
            : m_num_nodes(num_nodes), m_size(size), m_iter_max(iter_max), m_range_min_list(range_min_list), m_range_max_list(range_max_list)
    {
        m_name = "DOA";
        m_P = 0.5;
        m_Q = 0.7;
        m_na_ini = 2;
        m_na_end = std::floor(m_size / 2.0);
    }

    void init();
    void run();
    void update_survival_rate();
    void update_position();
    std::vector<double> get_sumatory(int na, int id);
    int get_best_id();
    //void export_topology(const std::string& filename);
    std::vector<std::vector<bool>> export_topology();
protected:
    std::string m_name;
    double m_P, m_Q;
    int m_na_ini, m_na_end;
    int m_num_nodes, m_size, m_iter_max;
    std::vector<double> m_range_min_list, m_range_max_list;
    std::vector<DOA_Unit> m_unit_list;
    std::vector<double> m_position_best;

private:
    double cal_fitfunction(const std::vector<double>& position, const std::vector<std::vector<bool>>& connections);
    std::vector<double> get_out_bound_value(const std::vector<double>& position);
    std::vector<double> levy(int d);
};

void DOA_Base::init()
{
    std::random_device rd;
    std::mt19937 gen(42);
    for (int i = 0; i < m_size; ++i)
    {
        DOA_Unit unit;
        unit.position.resize(m_num_nodes * 2); // 二维位置
        unit.connections.resize(m_num_nodes, std::vector<bool>(m_num_nodes, false));

        // 随机生成节点位置
        for (int j = 0; j < m_num_nodes * 2; ++j)
        {
            std::uniform_real_distribution<> dis(m_range_min_list[j % 2], m_range_max_list[j % 2]);
            unit.position[j] = dis(gen);
        }

        // 随机生成连接
        for (int j = 0; j < m_num_nodes; ++j)
        {
            for (int k = j + 1; k < m_num_nodes; ++k)
            {
                if (std::uniform_real_distribution<>(0, 1)(gen) < 0.5)
                {
                    unit.connections[j][k] = unit.connections[k][j] = true;
                }
            }
        }

        unit.value = cal_fitfunction(unit.position, unit.connections);
        unit.survival_rate = 0.0; // 初始化为0
        m_unit_list.push_back(unit);
    }
    m_position_best = m_unit_list[get_best_id()].position;
}

void DOA_Base::run()
{
    for (int iter = 0; iter < m_iter_max; ++iter)
    {
        update_survival_rate();
        update_position();
        std::cout << "Iteration " << iter << ": Best Value = " << m_unit_list[get_best_id()].value << std::endl;
    }
}

void DOA_Base::update_survival_rate()
{
    std::vector<double> values;
    for (const auto& unit : m_unit_list)
    {
        values.push_back(unit.value);
    }
    std::sort(values.begin(), values.end(), std::greater<double>());

    for (auto& unit : m_unit_list)
    {
        double survival_rate = std::abs(values[0] - unit.value) / std::abs(values[0] - values.back()) + std::numeric_limits<double>::min();
        unit.survival_rate = survival_rate;
    }
}

void DOA_Base::update_position()
{
    std::random_device rd;
    std::mt19937 gen(42);
    std::uniform_int_distribution<> dis_na(m_na_ini, m_na_end);
    std::uniform_real_distribution<> dis(-2, 2);
    std::uniform_real_distribution<> dis1(-1, 1);
    std::uniform_real_distribution<> dis_scavenging(-1, 1);
    std::uniform_real_distribution<> dis_binary(0, 1);

    for (int i = 0; i < m_size; ++i)
    {
        int na = dis_na(gen);
        if (std::uniform_real_distribution<>(0, 1)(gen) < m_P)
        {
            if (std::uniform_real_distribution<>(0, 1)(gen) < m_Q)
            {
                auto pos_sumatory = get_sumatory(na, i);
                std::vector<double> new_pos;
                for (int j = 0; j < m_num_nodes * 2; ++j)
                {
                    new_pos.push_back(pos_sumatory[j] * dis(gen) - m_position_best[j]);
                }
                new_pos = get_out_bound_value(new_pos);
                std::vector<std::vector<bool>> new_connections = m_unit_list[i].connections; // 保留原连接
                double new_value = cal_fitfunction(new_pos, new_connections);
                if (new_value > m_unit_list[i].value)
                {
                    m_unit_list[i].position = new_pos;
                    m_unit_list[i].value = new_value;
                }
            }
            else
            {
                // 这里简化了逻辑，实际上应包括连接的更新
                int r_id = std::uniform_int_distribution<>(0, m_size - 1)(gen);
                double beta1 = dis(gen);
                double beta2 = dis1(gen);
                std::vector<double> new_pos;
                for (int j = 0; j < m_num_nodes * 2; ++j)
                {
                    new_pos.push_back(m_position_best[j] + beta1 * std::exp(beta2) * (m_unit_list[r_id].position[j] - m_unit_list[i].position[j]));
                }
                new_pos = get_out_bound_value(new_pos);
                std::vector<std::vector<bool>> new_connections = m_unit_list[i].connections; // 保留原连接
                double new_value = cal_fitfunction(new_pos, new_connections);
                if (new_value > m_unit_list[i].value)
                {
                    m_unit_list[i].position = new_pos;
                    m_unit_list[i].value = new_value;
                }
            }
        }
        else
        {
            // 食腐行为
            int r_id = std::uniform_int_distribution<>(0, m_size - 1)(gen);
            double m = dis_scavenging(gen);
            double n = dis_binary(gen) <= 0.5 ? 0 : 1;
            std::vector<double> new_pos;
            for (int j = 0; j < m_num_nodes * 2; ++j)
            {
                new_pos.push_back(0.5 * (std::exp(m) * m_unit_list[r_id].position[j] - std::pow(-1, n) * m_unit_list[i].position[j]));
            }
            new_pos = get_out_bound_value(new_pos);
            std::vector<std::vector<bool>> new_connections = m_unit_list[i].connections; // 保留原连接
            double new_value = cal_fitfunction(new_pos, new_connections);
            if (new_value > m_unit_list[i].value)
            {
                m_unit_list[i].position = new_pos;
                m_unit_list[i].value = new_value;
            }
        }

        // 对于存活率低于0.3的个体，使用特殊策略更新位置
        if (m_unit_list[i].survival_rate < 0.3)
        {
            int r_id1 = std::uniform_int_distribution<>(0, m_size - 1)(gen);
            int r_id2 = std::uniform_int_distribution<>(0, m_size - 1)(gen);
            double binary = (dis_binary(gen) < 0.5) ? -1.0 : 1.0;
            std::vector<double> new_pos;
            for (int j = 0; j < m_num_nodes * 2; ++j)
            {
                new_pos.push_back(m_position_best[j] + (m_unit_list[r_id1].position[j] + binary * m_unit_list[r_id2].position[j]) / 2);
            }
            new_pos = get_out_bound_value(new_pos);
            std::vector<std::vector<bool>> new_connections = m_unit_list[i].connections; // 保留原连接
            double new_value = cal_fitfunction(new_pos, new_connections);
            if (new_value > m_unit_list[i].value)
            {
                m_unit_list[i].position = new_pos;
                m_unit_list[i].value = new_value;
            }
        }
    }
}

std::vector<double> DOA_Base::get_sumatory(int na, int id)
{
    std::vector<int> A(m_size);
    std::iota(A.begin(), A.end(), 0);
    std::mt19937 gen(42); // 添加这一行
    std::random_shuffle(A.begin(), A.end(), [&gen](int i) { return gen() % (i + 1); });
    std::vector<int> r_ids;
    for (int i = 0; i < na + 1; ++i)
    {
        if (A[i] != id)
            r_ids.push_back(A[i]);
    }

    std::vector<double> pos_sumatory(m_num_nodes * 2, 0.0);
    for (auto rid : r_ids)
    {
        for (int i = 0; i < m_num_nodes * 2; ++i)
        {
            pos_sumatory[i] += m_unit_list[id].position[i] - m_unit_list[rid].position[i];
        }
    }
    for (auto& p : pos_sumatory)
    {
        p /= r_ids.size();
    }
    return pos_sumatory;
}

int DOA_Base::get_best_id()
{
    return std::max_element(m_unit_list.begin(), m_unit_list.end(),
                            [](const DOA_Unit& a, const DOA_Unit& b) { return a.value < b.value; }) - m_unit_list.begin();
}

// 适应度函数
double DOA_Base::cal_fitfunction(const std::vector<double>& position, const std::vector<std::vector<bool>>& connections)
{
    // 这里简化了适应度计算，实际应用中需要根据具体的网络拓扑需求来实现
    double fitness = 0.0;
    for (int i = 0; i < m_num_nodes; ++i)
    {
        for (int j = i + 1; j < m_num_nodes; ++j)
        {
            if (connections[i][j])
            {
                double distance = std::sqrt(std::pow(position[i * 2] - position[j * 2], 2) + std::pow(position[i * 2 + 1] - position[j * 2 + 1], 2));
                fitness += 1.0 / distance; // 假设更短的距离带来更高的适应度
            }
        }
    }
    return fitness;
}

std::vector<double> DOA_Base::get_out_bound_value(const std::vector<double>& position)
{
    std::vector<double> new_pos = position;
    for (int i = 0; i < m_num_nodes * 2; ++i)
    {
        if (new_pos[i] < m_range_min_list[i % 2] || new_pos[i] > m_range_max_list[i % 2])
        {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> dis(m_range_min_list[i % 2], m_range_max_list[i % 2]);
            new_pos[i] = dis(gen);
        }
    }
    return new_pos;
}

// Levy函数的实现
std::vector<double> DOA_Base::levy(int d)
{
    double beta = 1.5;
    double sigma = std::pow((std::tgamma(1 + beta) * std::sin(M_PI * beta / 2) / (std::tgamma((1 + beta) / 2) * beta * std::pow(2, (beta - 1) / 2))), 1 / beta);
    std::random_device rd;
    std::mt19937 gen(42);
    std::normal_distribution<> normal(0, sigma);

    std::vector<double> step(d);
    for (int i = 0; i < d; ++i)
    {
        double u = normal(gen);
        double v = std::normal_distribution<>(0, 1)(gen);
        step[i] = u / std::pow(std::abs(v), 1 / beta);
    }
    return step;
}



std::vector<Edge> topologyToEdges(const std::vector<std::vector<bool>>& topology) {
    std::vector<Edge> edges;

    for (size_t i = 0; i < topology.size(); ++i) {
        for (size_t j = i + 1; j < topology[i].size(); ++j) {
            if (topology[i][j]) {
                bool duplicate = false;
                for (const auto& e : edges) {
                    if ((e.src == static_cast<int>(i) && e.dest == static_cast<int>(j)) ||
                        (e.src == static_cast<int>(j) && e.dest == static_cast<int>(i))) {
                        duplicate = true;
                        break;
                    }
                }
                if (!duplicate) {
                    edges.push_back({static_cast<int>(i), static_cast<int>(j)});
                }
            }
        }
    }
    return edges;
}


std::vector<std::vector<bool>> DOA_Base::export_topology() {
    std::vector<std::vector<bool>> topology(m_num_nodes, std::vector<bool>(m_num_nodes, false));
    int best_id = get_best_id();
    for(int i = 0; i < m_num_nodes; i++) {
        for(int j = i + 1; j < m_num_nodes; j++) { // 只考虑上三角矩阵
            topology[i][j] = topology[j][i] = m_unit_list[best_id].connections[i][j];
        }
    }
    return topology;
}


int main()
{
    int num_nodes = 25;
    int size = 50;
    int iter_max = 200;
    std::vector<double> range_min_list = {-20, -20}; // 节点位置范围
    std::vector<double> range_max_list = {20, 20};


    DOA_Base doa(num_nodes, size, iter_max, range_min_list, range_max_list);
    doa.init();
    doa.run();
    std::vector<Edge> edges = topologyToEdges(doa.export_topology());

      // 创建节点
    NodeContainer nodes;
    nodes.Create(num_nodes);

//    // 安装互联网堆栈（IPv4）
    InternetStackHelper internet;
    internet.Install(nodes);

    // 点对点助手
    PointToPointHelper pointToPoint;
    pointToPoint.SetDeviceAttribute("DataRate", StringValue("5Mbps"));
    pointToPoint.SetChannelAttribute("Delay", StringValue("2ms"));
//
//   std::map<int, Ipv4InterfaceContainer> interfaces;
//int subnetCounter = 0;
// int MAX_SUBNETS = 1024; // 根据你的网络设计调整这个值
//for (auto edge : edges) {
//    NodeContainer link = NodeContainer(nodes.Get(edge.src), nodes.Get(edge.dest));
//    NetDeviceContainer devices = pointToPoint.Install(link);
//    Ipv4AddressHelper address;
//    if (subnetCounter >= MAX_SUBNETS) {
//        NS_FATAL_ERROR("Exceeded maximum number of subnets.");
//    }
//    std::string baseAddress = "10.1." + std::to_string(subnetCounter++) + ".0";
//    address.SetBase(baseAddress.c_str(), "255.255.255.0");
//    interfaces[edge.src] = address.Assign(devices);
//}
std::map<int, Ipv4InterfaceContainer> interfaces;
std::set<int> allocatedNodes; // 记录已经分配了IP地址的节点
int subnetCounter = 0;
int MAX_SUBNETS = 1024;

for (auto edge : edges) {
    if (allocatedNodes.find(edge.src) == allocatedNodes.end()) {
        NodeContainer link = NodeContainer(nodes.Get(edge.src), nodes.Get(edge.dest));
        NetDeviceContainer devices = pointToPoint.Install(link);
        Ipv4AddressHelper address;
        if (subnetCounter >= MAX_SUBNETS) {
            NS_FATAL_ERROR("Exceeded maximum number of subnets.");
        }
        std::string baseAddress = "10.1." + std::to_string(subnetCounter++) + ".0";
        address.SetBase(baseAddress.c_str(), "255.255.255.0");
        interfaces[edge.src] = address.Assign(devices);
        allocatedNodes.insert(edge.src);
    }
    // 对于目的节点的处理逻辑类似
}


  Ipv4GlobalRoutingHelper::PopulateRoutingTables();

    // 设置仿真时间
    Simulator::Stop(Seconds(10.0));

    // 为网络中的每条边安装应用
    uint16_t port = 9; // 监听端口
    for (auto edge : edges) {
        // 服务器端（使用PacketSink）
        PacketSinkHelper sinkHelper("ns3::UdpSocketFactory", InetSocketAddress(interfaces[edge.src].GetAddress(1), port));
        ApplicationContainer sinkApps = sinkHelper.Install(nodes.Get(edge.dest));
        sinkApps.Start(Seconds(1.0));
        sinkApps.Stop(Seconds(10.0));

        // 客户端（使用OnOffApplication）
        OnOffHelper onoff("ns3::UdpSocketFactory", Address(InetSocketAddress(interfaces[edge.src].GetAddress(1), port)));
        onoff.SetConstantRate(DataRate("5Mbps"));
        ApplicationContainer apps = onoff.Install(nodes.Get(edge.src));
        apps.Start(Seconds(1.1));
        apps.Stop(Seconds(10.0));
    }

      FlowMonitorHelper flowmon;
    Ptr<FlowMonitor> monitor = flowmon.InstallAll();

    // 开始计时
    auto start = std::chrono::high_resolution_clock::now();

    // 运行仿真
    Simulator::Run();
    Simulator::Destroy();

    // 结束计时
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    std::cout << "Total simulation time: " << elapsed.count() << " ms" << std::endl;


     monitor->CheckForLostPackets();
    Ptr<Ipv4FlowClassifier> classifier = DynamicCast<Ipv4FlowClassifier>(flowmon.GetClassifier());
    std::map<FlowId, FlowMonitor::FlowStats> stats = monitor->GetFlowStats();

//    double avgThroughPut = 0.0;
//    double avgDelay = 0.0;
//    int sss = 0;
//    for (std::map<FlowId, FlowMonitor::FlowStats>::const_iterator i = stats.begin(); i != stats.end(); ++i) {
//    sss++;
//        Ipv4FlowClassifier::FiveTuple t = classifier->FindFlow(i->first);
//        std::cout << "Flow " << i->first << " (" << t.sourceAddress << " -> " << t.destinationAddress << ")\n";
//        std::cout << "  Tx Bytes: " << i->second.txBytes << "\n";
//        std::cout << "  Rx Bytes: " << i->second.rxBytes << "\n";
//        std::cout << "  Throughput: " << i->second.rxBytes * 8.0 / (10.0 - 1.1) / 1e6 << " Mbps\n";
//        avgThroughPut += i->second.rxBytes * 8.0 / (10.0 - 1.1) / 1e6;
//        std::cout << "  Delay: " << i->second.delaySum.GetSeconds() / i->second.rxPackets << " seconds\n";
//        avgDelay += i->second.delaySum.GetSeconds() / i->second.rxPackets;
//    }
//    std::cout<<"sss: "<<sss <<"\n";
//    std::cout<< "avgThroughPut: " << avgThroughPut / sss << " Mbps\n";
//    std::cout<< "avgDelay: " << avgDelay / sss << " seconds\n";
//
//    return 0;
 // 收集每个节点的流数据
  double totalThroughput = 0.0;
    double totalDelay = 0.0;
    int totalFlows = 0;

    for (std::map<FlowId, FlowMonitor::FlowStats>::const_iterator i = stats.begin(); i != stats.end(); ++i) {
        totalThroughput += i->second.rxBytes * 8.0 / (10.0 - 1.1) / 1e6; // 累加吞吐量
        totalDelay += i->second.delaySum.GetSeconds() / i->second.rxPackets; // 累加时延
        ++totalFlows; // 流的数量
    }

    // 计算总平均吞吐量和总平均时延
    if (totalFlows > 0) {
        double avgThroughput = totalThroughput / num_nodes;
        double avgDelay = totalDelay / num_nodes;

        std::cout << "Total Average Throughput: " << avgThroughput << " Mbps" << std::endl;
        std::cout << "Total Average Delay: " << avgDelay << " seconds" << std::endl;
    } else {
        std::cout << "No flows detected." << std::endl;
    }

    return 0;
}
