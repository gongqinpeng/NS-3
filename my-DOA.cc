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
#include <cmath>
#include <limits>

using namespace ns3;

struct Edge {
    int src, dest, weight;
};

struct Dingo {
    std::vector<std::vector<bool>> topology; // 邻接矩阵表示的拓扑
    std::vector<double> velocity; // 速度向量
    double fitness; // 适应度值
    std::vector<std::vector<bool>> pBestTopology; // 个人最佳拓扑
    double pBestFitness; // 个人最佳适应度值
};

// 计算拓扑的适应度，这里假设我们要最小化拓扑中的总权重

double calculateFitness(const std::vector<std::vector<bool>>& topology) {
    double totalThroughput = 0.0;
    double totalDelay = 0.0;
    std::random_device rd;
    std::mt19937 gen(42);
    std::uniform_real_distribution<> disThroughput(0.1, 5.0); // 吞吐量范围
    std::uniform_real_distribution<> disDelay(1.0, 5.0); // 延迟范围

    for (size_t i = 0; i < topology.size(); ++i) {
        for (size_t j = i + 1; j < topology[i].size(); ++j) {
            if (topology[i][j]) {
                double throughput = disThroughput(gen);
                double delay = disDelay(gen);

                totalThroughput += throughput;
                totalDelay += delay;
            }
        }
    }

    // 适应度为吞吐量的加权和与延迟的倒数
    return (totalThroughput / (topology.size() * (topology.size() - 1))) / totalDelay;

}

// 初始化野狗群
std::vector<Dingo> initializePack(int numNodes, int numDingos) {
    std::vector<Dingo> pack(numDingos);
    std::random_device rd;
    std::mt19937 gen(42); // 使用固定种子以保证结果的一致性
    std::uniform_int_distribution<> dis(0, numNodes - 1);

    for (auto& dingo : pack) {
        // 初始化野狗的拓扑
        dingo.topology = std::vector<std::vector<bool>>(numNodes, std::vector<bool>(numNodes, false));
        // 随机连接一些节点
        for (int i = 0; i < numNodes; ++i) {
            int j = dis(gen);
            if (i != j) {
                dingo.topology[i][j] = dingo.topology[j][i] = true;
            }
        }
        // 初始化速度和适应度
        dingo.velocity = std::vector<double>(numNodes * numNodes, 5.0);
        dingo.fitness = calculateFitness(dingo.topology);
        dingo.pBestTopology = dingo.topology;
        dingo.pBestFitness = dingo.fitness;
    }
    return pack;
}



// DOA优化算法
std::vector<std::vector<bool>> doaOptimization(int numNodes, int numDingos, int maxIterations, double P, double Q) {
    auto pack = initializePack(numNodes, numDingos);
    std::vector<std::vector<bool>> gBestTopology = pack[0].topology;
    double gBestFitness = pack[0].fitness;

    std::random_device rd;
    std::mt19937 gen(42); // 使用固定种子以保证结果的一致性
    std::uniform_real_distribution<> dis(0.0, 1.0);

    for (int iter = 0; iter < maxIterations; ++iter) {
        for (auto& dingo : pack) {
            double r = dis(gen);
            if (r < P) { // 如果随机数小于P
                if (dis(gen) < Q) { // 进入到Q
                    // 围攻行为：尝试向最佳解靠近
                    for (int i = 0; i < numNodes; ++i) {
                        for (int j = 0; j < numNodes; ++j) {
                            if (i != j) {
                                double r1 = dis(gen);
                                double r2 = dis(gen);
                                double v = dingo.velocity[i * numNodes + j];
                                double cognitive = 2.0 * r1 * (dingo.pBestTopology[i][j] - dingo.topology[i][j]);
                                double social = 2.0 * r2 * (gBestTopology[i][j] - dingo.topology[i][j]);
                                dingo.velocity[i * numNodes + j] = v + cognitive + social;

                                // 更新拓扑
                                dingo.topology[i][j] = dingo.topology[j][i] = (dis(gen) < std::exp(-dingo.velocity[i * numNodes + j])) ? true : false;
                            }
                        }
                    }
                } else { // 大于Q
                    // 捕食行为：可能进行更激进的探索
                    // 这里可以是更大的移动或随机搜索
                    for (int i = 0; i < numNodes; ++i) {
                        for (int j = 0; j < numNodes; ++j) {
                            if (i != j && dis(gen) < 0.1) { // 10%的概率改变连接
                                dingo.topology[i][j] = dingo.topology[j][i] = !dingo.topology[i][j];
                            }
                        }
                    }
                }
            } else { // 大于P
                // 食腐行为：采用其他野狗的最佳解或随机改变拓扑
                if (dis(gen) < 0.5) { // 50%的概率采用其他野狗的最佳解
                    auto otherDingo = pack[dis(gen) * numDingos]; // 随机选择另一只野狗
                    dingo.topology = otherDingo.pBestTopology;
                } else { // 50%的概率随机改变拓扑
                    for (int i = 0; i < numNodes; ++i) {
                        for (int j = 0; j < numNodes; ++j) {
                            if (i != j && dis(gen) < 0.1) { // 10%的概率改变连接
                                dingo.topology[i][j] = dingo.topology[j][i] = !dingo.topology[i][j];
                            }
                        }
                    }
                }
            }

            // 计算新适应度
            dingo.fitness = calculateFitness(dingo.topology);

            // 更新个人最佳
            if (dingo.fitness > dingo.pBestFitness) {
                dingo.pBestTopology = dingo.topology;
                dingo.pBestFitness = dingo.fitness;
            }

            // 更新全局最佳
            if (dingo.fitness > gBestFitness) {
                gBestTopology = dingo.topology;
                gBestFitness = dingo.fitness;
            }

            // 生存：确保适应度不低于某一阈值
            if (dingo.fitness < 0.01) { // 假设0.01是一个低适应度阈值
                dingo = initializePack(numNodes, 1)[0]; // 重新初始化野狗
            }
        }
    }

    return gBestTopology;
}
// 将邻接矩阵转换为边的列表
std::vector<Edge> topologyToEdges(const std::vector<std::vector<bool>>& topology) {
    std::vector<Edge> edges;
    std::random_device rd;
    std::mt19937 gen(42); // 使用固定种子以保证结果的一致性
    std::uniform_int_distribution<> dis(1, 10); // 权重范围为1到10

    for (size_t i = 0; i < topology.size(); ++i) {
        for (size_t j = i + 1; j < topology[i].size(); ++j) {
            if (topology[i][j]) {
                edges.push_back({static_cast<int>(i), static_cast<int>(j), dis(gen)});
            }
        }
    }
    return edges;
}

void exportTopologyToJson(const std::vector<std::vector<bool>>& topology, const std::string& filename) {
    std::ofstream jsonFile(filename);
    jsonFile << "[" << std::endl;
    for (size_t i = 0; i < topology.size(); ++i) {
        jsonFile << "    [";
        for (size_t j = 0; j < topology[i].size(); ++j) {
            jsonFile << (topology[i][j] ? "1" : "0");
            if (j < topology[i].size() - 1) jsonFile << ",";
        }
        jsonFile << "]";
        if (i < topology.size() - 1) jsonFile << ",";
        jsonFile << std::endl;
    }
    jsonFile << "]" << std::endl;
}

int main(int argc, char *argv[]) {
    CommandLine cmd;
    int numNodes = 80; // 默认节点数量
    int numEdges = 190; // 默认边的数量
    int numDingos = 50; // DOA的野狗数量
    int maxIterations = 100; // DOA的最大迭代次数
double P = 0.5; // 围攻和捕食的概率
    double Q = 0.7; // 围攻的概率
    // 命令行参数解析
    cmd.AddValue("numNodes", "Number of nodes", numNodes);
    cmd.AddValue("numEdges", "Number of edges", numEdges);
    cmd.AddValue("numDingos", "Number of dingos in DOA", numDingos);
    cmd.AddValue("maxIterations", "Maximum iterations for DOA", maxIterations);
    cmd.Parse(argc, argv);

    // 生成网络拓扑
    std::vector<std::vector<bool>> topology = doaOptimization(numNodes, numDingos, maxIterations,P,Q);
    std::vector<Edge> edges = topologyToEdges(topology);

    exportTopologyToJson(topology, "/usr/womensan/Desktop/ns-3-dev/scratch/DOA-80-topology.json");

    // 创建节点
//    NodeContainer nodes;
//    nodes.Create(numNodes);
//
////    // 安装互联网堆栈（IPv4）
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
//    for (auto edge : edges) {
//        NodeContainer link = NodeContainer(nodes.Get(edge.src), nodes.Get(edge.dest));
//        NetDeviceContainer devices = pointToPoint.Install(link);
//        Ipv4AddressHelper address;
//        std::string baseAddress = "10.1." + std::to_string(subnetCounter++) + ".0"; // 使用唯一子网号
//        address.SetBase(baseAddress.c_str(), "255.255.255.0");
//        interfaces[edge.src] = address.Assign(devices);
//    }
//
//    // 配置路由
//  Ipv4GlobalRoutingHelper::PopulateRoutingTables();
//
    // 设置仿真时间
    Simulator::Stop(Seconds(10.0));
//
//    // 为网络中的每条边安装应用
//    uint16_t port = 9; // 监听端口
//    for (auto edge : edges) {
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
//
//    // 配置FlowMonitor
//    FlowMonitorHelper flowmon;
//    Ptr<FlowMonitor> monitor = flowmon.InstallAll();
//
    // 开始计时
    auto start = std::chrono::high_resolution_clock::now();

    // 运行仿真
    Simulator::Run();
    Simulator::Destroy();

    // 结束计时
    auto end = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double, std::milli> elapsed = end - start;
//    std::cout << "Total simulation time: " << elapsed.count() << " ms" << std::endl;
//
////    // 收集流量监控数据
//    monitor->CheckForLostPackets();
//    Ptr<Ipv4FlowClassifier> classifier = DynamicCast<Ipv4FlowClassifier>(flowmon.GetClassifier());
//    std::map<FlowId, FlowMonitor::FlowStats> stats = monitor->GetFlowStats();
//
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

    return 0;
}
