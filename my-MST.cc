#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/point-to-point-module.h"
#include "ns3/applications-module.h"
#include "ns3/ipv4-flow-classifier.h"
#include "ns3/flow-monitor-module.h"
#include <vector>
#include <map>
#include <iostream>
#include <random>
#include <algorithm>
#include <chrono>

using namespace ns3;

struct Edge {
    int src, dest, weight;
};

// 生成随机图，并确保边不重复
std::vector<Edge> generateRandomGraph(int V, int E) {

    std::vector<Edge> graph;
    std::random_device rd;
    std::mt19937 gen(42);
    std::uniform_int_distribution<> dis(1, 10); // 权重范围为1到10

    for (int i = 0; i < E && graph.size() < V*(V-1)/2; ++i) {
        int src = dis(gen) % V;
        int dest = (src + dis(gen) % (V - 1) + 1) % V; // 确保不连接到自身
        int weight = dis(gen);

        // 检查是否已经存在这条边，如果存在则更新权重
        auto it = std::find_if(graph.begin(), graph.end(), [&](const Edge& e) {
            return (e.src == src && e.dest == dest) || (e.src == dest && e.dest == src);
        });
        if (it == graph.end()) {
            graph.push_back({src, dest, weight});
        } else {
            // 如果边已经存在，选择权重较小的边
            if (weight < it->weight) {
                *it = {src, dest, weight};
            }
        }
    }
    return graph;
}

// Prim MST算法
std::vector<Edge> primMST(std::vector<Edge>& edges, int V) {
    std::vector<int> parent(V, -1);
    std::vector<int> key(V, INT_MAX);
    std::vector<bool> mstSet(V, false);

    key[0] = 0;
    parent[0] = -1;

    for (int count = 0; count < V - 1; count++) {
        int u = -1;
        for (int v = 0; v < V; v++)
            if (!mstSet[v] && (u == -1 || key[v] < key[u]))
                u = v;

        mstSet[u] = true;

        for (auto& edge : edges) {
            if (edge.src == u || edge.dest == u) {
                int v = (edge.src == u) ? edge.dest : edge.src;
                if(v >= V) continue;
                int weight = edge.weight;

                if (!mstSet[v] && key[v] > weight) {
                    key[v] = weight;
                    parent[v] = u;
                }
            }
        }
    }

    std::vector<Edge> result;
    for (int i = 1; i < V; ++i)
       if (parent[i] != -1) // 确保parent[i]有效
            result.push_back({parent[i], i, key[i]});

    return result;
}

int main(int argc, char *argv[]) {
    CommandLine cmd;
    int numNodes = 100; // 默认节点数量
    int numEdges = 4950; // 默认边的数量

    // 命令行参数解析
    cmd.AddValue("numNodes", "Number of nodes", numNodes);
    cmd.AddValue("numEdges", "Number of edges", numEdges);
    cmd.Parse(argc, argv);
      auto start = std::chrono::high_resolution_clock::now();
    // 生成随机图
    std::vector<Edge> graph = generateRandomGraph(numNodes, numEdges);

    // 生成MST
    std::vector<Edge> mst = primMST(graph, numNodes);
      auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    std::cout << "Total simulation time: " << elapsed.count() << " ms" << std::endl;
//    // 创建节点
//    NodeContainer nodes;
//    nodes.Create(numNodes);
//
//    // 安装互联网堆栈（IPv4）
//    InternetStackHelper internet;
//    internet.Install(nodes);
//
//    // 点对点助手
//    PointToPointHelper pointToPoint;
//    pointToPoint.SetDeviceAttribute("DataRate", StringValue("5Mbps"));
//    pointToPoint.SetChannelAttribute("Delay", StringValue("2ms"));
//
//    // 安装设备并分配IP地址
//    std::map<int, Ipv4InterfaceContainer> interfaces;
//    int subnetCounter = 0;
//    for (auto edge : mst) {
//        NodeContainer link = NodeContainer(nodes.Get(edge.src), nodes.Get(edge.dest));
//        NetDeviceContainer devices = pointToPoint.Install(link);
//        Ipv4AddressHelper address;
//        std::string baseAddress = "10.1." + std::to_string(subnetCounter++) + ".0"; // 使用唯一子网号
//        address.SetBase(baseAddress.c_str(), "255.255.255.0");
//        interfaces[edge.src] = address.Assign(devices);
//    }
//
//    // 配置路由
//    Ipv4GlobalRoutingHelper::PopulateRoutingTables();
//
//    // 设置仿真时间
//    Simulator::Stop(Seconds(10.0));
//
//    // 为MST中的每条边安装应用
//    uint16_t port = 9; // 监听端口
//    for (auto edge : mst) {
//        // 服务器端（使用PacketSink）
//        PacketSinkHelper sinkHelper("ns3::UdpSocketFactory", InetSocketAddress(interfaces[edge.src].GetAddress(1), port));
//        ApplicationContainer sinkApps = sinkHelper.Install(nodes.Get(edge.dest));
//        sinkApps.Start(Seconds(1.0));
//        sinkApps.Stop(Seconds(10.0));
//
//        // 客户端（使用OnOffApplication）
//        OnOffHelper onoff("ns3::UdpSocketFactory", Address(InetSocketAddress(interfaces[edge.src].GetAddress(1), port)));
//        onoff.SetConstantRate(DataRate("5Mbps"));
//        ApplicationContainer apps = onoff.Install(nodes.Get(edge.src));
//        apps.Start(Seconds(1.1));
//        apps.Stop(Seconds(10.0));
//    }
////
////     //配置FlowMonitor
//    FlowMonitorHelper flowmon;
//    Ptr<FlowMonitor> monitor = flowmon.InstallAll();
//
//    // 开始计时
//    auto start = std::chrono::high_resolution_clock::now();
//
//    // 运行仿真
//    Simulator::Run();
//    Simulator::Destroy();
//
//    // 结束计时
//    auto end = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double, std::milli> elapsed = end - start;
//std::cout << "Total simulation time: " << elapsed.count() << " ms" << std::endl;
//    // 收集流量监控数据
//    monitor->CheckForLostPackets();
//    Ptr<Ipv4FlowClassifier> classifier = DynamicCast<Ipv4FlowClassifier>(flowmon.GetClassifier());
//    std::map<FlowId, FlowMonitor::FlowStats> stats = monitor->GetFlowStats();
////
////    double avgThroughPut = 0.0;
////    double avgDelay = 0.0;
////    int sss = 0;
////    for (std::map<FlowId, FlowMonitor::FlowStats>::const_iterator i = stats.begin(); i != stats.end(); ++i) {
////    sss++;
////        Ipv4FlowClassifier::FiveTuple t = classifier->FindFlow(i->first);
////        std::cout << "Flow " << i->first << " (" << t.sourceAddress << " -> " << t.destinationAddress << ")\n";
////        std::cout << "  Tx Bytes: " << i->second.txBytes << "\n";
////        std::cout << "  Rx Bytes: " << i->second.rxBytes << "\n";
////        std::cout << "  Throughput: " << i->second.rxBytes * 8.0 / (10.0 - 1.1) / 1e6 << " Mbps\n";
////        avgThroughPut += i->second.rxBytes * 8.0 / (10.0 - 1.1) / 1e6;
////        std::cout << "  Delay: " << i->second.delaySum.GetSeconds() / i->second.rxPackets << " seconds\n";
////        avgDelay += i->second.delaySum.GetSeconds() / i->second.rxPackets;
////    }
//// std::cout<< "sss: " << sss <<"\n";
////    std::cout<< "avgThroughPut: " << avgThroughPut / numNodes << " Mbps\n";
////    std::cout<< "avgDelay: " << avgDelay / numNodes << " seconds\n";
//
//  double totalThroughput = 0.0;
//    double totalDelay = 0.0;
//    int totalFlows = 0;
//
//    for (std::map<FlowId, FlowMonitor::FlowStats>::const_iterator i = stats.begin(); i != stats.end(); ++i) {
//        totalThroughput += i->second.rxBytes * 8.0 / (10.0 - 1.1) / 1e6; // 累加吞吐量
//        totalDelay += i->second.delaySum.GetSeconds() / i->second.rxPackets; // 累加时延
//        ++totalFlows; // 流的数量
//    }
//
//    // 计算总平均吞吐量和总平均时延
//    if (totalFlows > 0) {
//        double avgThroughput = totalThroughput / numNodes;
//        double avgDelay = totalDelay / numNodes;
//
//        std::cout << "Total Average Throughput: " << avgThroughput << " Mbps" << std::endl;
//        std::cout << "Total Average Delay: " << avgDelay << " seconds" << std::endl;
//    } else {
//        std::cout << "No flows detected." << std::endl;
//    }

    return 0;
}

