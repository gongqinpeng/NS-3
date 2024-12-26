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

struct Particle {
    std::vector<std::vector<bool>> topology; // 邻接矩阵表示的拓扑
    std::vector<double> velocity; // 速度向量
    double fitness; // 适应度值
    std::vector<std::vector<bool>> pBestTopology; // 个人最佳拓扑
    double pBestFitness; // 个人最佳适应度值
};


// 计算拓扑的适应度，这里假设我们要最小化拓扑中的总权重
double calculateFitness(const std::vector<std::vector<bool>>& topology) {
    double totalWeight = 0.0;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(1, 10); // 权重范围为1到10

    for (size_t i = 0; i < topology.size(); ++i) {
        for (size_t j = i + 1; j < topology[i].size(); ++j) {
            if (topology[i][j]) {
                totalWeight += dis(gen); // 随机权重
            }
        }
    }
    return 1.0 / totalWeight; // 适应度为权重的倒数
}
// 初始化粒子群
std::vector<Particle> initializeSwarm(int numNodes, int numParticles) {
    std::vector<Particle> swarm(numParticles);
    std::random_device rd;
    std::mt19937 gen(42); // 使用固定种子以保证结果的一致性
    std::uniform_int_distribution<> dis(0, numNodes - 1);

    for (auto& particle : swarm) {
        // 初始化粒子的拓扑
        particle.topology = std::vector<std::vector<bool>>(numNodes, std::vector<bool>(numNodes, false));
        // 随机连接一些节点
        for (int i = 0; i < numNodes; ++i) {
            int j = dis(gen);
            if (i != j) {
                particle.topology[i][j] = particle.topology[j][i] = true;
            }
        }
        // 初始化速度和适应度
        particle.velocity = std::vector<double>(numNodes * numNodes, 0.0);
        particle.fitness = calculateFitness(particle.topology);
        particle.pBestTopology = particle.topology;
        particle.pBestFitness = particle.fitness;
    }
    return swarm;
}



// PSO优化
std::vector<std::vector<bool>> optimizeTopology(int numNodes, int numParticles, int maxIterations) {
    auto swarm = initializeSwarm(numNodes, numParticles);
    std::vector<std::vector<bool>> gBestTopology = swarm[0].topology;
    double gBestFitness = swarm[0].fitness;

    std::random_device rd;
    std::mt19937 gen(42); // 使用固定种子以保证结果的一致性
    std::uniform_real_distribution<> dis(0.0, 1.0);

    for (int iter = 0; iter < maxIterations; ++iter) {
        for (auto& particle : swarm) {
            // 更新速度和位置
            for (int i = 0; i < numNodes; ++i) {
                for (int j = 0; j < numNodes; ++j) {
                    if (i != j) {
                        double r1 = dis(gen);
                        double r2 = dis(gen);
                        double v = particle.velocity[i * numNodes + j];
                        double cognitive = 2.0 * r1 * (particle.pBestTopology[i][j] - particle.topology[i][j]);
                        double social = 2.0 * r2 * (gBestTopology[i][j] - particle.topology[i][j]);
                        particle.velocity[i * numNodes + j] = v + cognitive + social;

                        // 更新拓扑
                        particle.topology[i][j] = particle.topology[j][i] = (dis(gen) < std::exp(-particle.velocity[i * numNodes + j])) ? true : false;
                    }
                }
            }

            // 计算新适应度
            particle.fitness = calculateFitness(particle.topology);

            // 更新个人最佳
            if (particle.fitness > particle.pBestFitness) {
                particle.pBestTopology = particle.topology;
                particle.pBestFitness = particle.fitness;
            }

            // 更新全局最佳
            if (particle.fitness > gBestFitness) {
                gBestTopology = particle.topology;
                gBestFitness = particle.fitness;
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

int main(int argc, char *argv[]) {
    CommandLine cmd;
    int numNodes = 10; // 默认节点数量
    int numEdges = 45; // 默认边的数量
    int numParticles = 50; // PSO的粒子数量
    int maxIterations = 100; // PSO的最大迭代次数

    // 命令行参数解析
    cmd.AddValue("numNodes", "Number of nodes", numNodes);
    cmd.AddValue("numEdges", "Number of edges", numEdges);
    cmd.AddValue("numParticles", "Number of particles in PSO", numParticles);
    cmd.AddValue("maxIterations", "Maximum iterations for PSO", maxIterations);
    cmd.Parse(argc, argv);

  auto start = std::chrono::high_resolution_clock::now();
    // 生成网络拓扑
    std::vector<std::vector<bool>> topology = optimizeTopology(numNodes, numParticles, maxIterations);
    std::vector<Edge> edges = topologyToEdges(topology);
 auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> elapsed = end - start;
    std::cout << "Total simulation time: " << elapsed.count() << " ms" << std::endl;
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
//     //配置路由
//  Ipv4GlobalRoutingHelper::PopulateRoutingTables();
//
//    // 设置仿真时间
//    Simulator::Stop(Seconds(10.0));
//
////    // 为网络中的每条边安装应用
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
//    std::cout << "Total simulation time: " << elapsed.count() << " ms" << std::endl;
////
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
////    std::cout<< "sss: " << sss <<"\n";
////    std::cout<< "avgThroughPut: " << avgThroughPut / numNodes << " Mbps\n";
////    std::cout<< "avgDelay: " << avgDelay / numNodes << " seconds\n";
// // 收集每个节点的流数据
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
