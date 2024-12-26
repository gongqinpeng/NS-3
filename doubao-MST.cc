#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/point-to-point-module.h"
#include "ns3/applications-module.h"
#include <iostream>
#include <vector>

// 用于生成随机的链路权重（这里简单示例，可按需替换更合理的权重生成方式）
double GenerateRandomWeight()
{
    return (double)(rand() % 10 + 1);  // 生成1到10之间的随机数作为权重
}

int main(int argc, char *argv[])
{
    // 启用NS3系统日志，便于查看运行过程中的相关信息
    ns3::LogComponentEnable("UdpEchoClientApplication", ns3::LOG_LEVEL_INFO);
    ns3::LogComponentEnable("UdpEchoServerApplication", ns3::LOG_LEVEL_INFO);

    // 创建节点容器，可根据需求调整节点数量
    ns3::NodeContainer nodes;
    nodes.Create(10);  // 创建8个节点

    // 存储链路信息的结构体，包含两个节点指针和链路权重
    struct LinkInfo {
        ns3::Ptr<ns3::Node> node1;
        ns3::Ptr<ns3::Node> node2;
        double weight;
    };
    std::vector<LinkInfo> allLinks;

    // 生成所有可能的链路（全连接形式），并给每条链路赋予随机权重
    for (size_t i = 0; i < nodes.GetN(); ++i) {
        for (size_t j = i + 1; j < nodes.GetN(); ++j) {
            LinkInfo link;
            link.node1 = nodes.Get(i);
            link.node2 = nodes.Get(j);
            link.weight = GenerateRandomWeight();
            allLinks.push_back(link);
        }
    }

    // 这里简单实现Prim算法来求最小生成树（实际可考虑优化及更完善的实现）
    std::vector<LinkInfo> mstLinks;
    std::vector<ns3::Ptr<ns3::Node>> visitedNodes;
    visitedNodes.push_back(nodes.Get(0));  // 从第一个节点开始

    while (visitedNodes.size() < nodes.GetN() && allLinks.size() > 0) {
        double minWeight = std::numeric_limits<double>::max();
        size_t minIndex = 0;
        for (size_t i = 0; i < allLinks.size(); ++i) {
            bool node1Visited = std::find(visitedNodes.begin(), visitedNodes.end(), allLinks[i].node1)!= visitedNodes.end();
            bool node2Visited = std::find(visitedNodes.begin(), visitedNodes.end(), allLinks[i].node2)!= visitedNodes.end();
            if ((node1Visited &&!node2Visited) || (!node1Visited && node2Visited)) {
                if (allLinks[i].weight < minWeight) {
                    minWeight = allLinks[i].weight;
                    minIndex = i;
                }
            }
        }
        mstLinks.push_back(allLinks[minIndex]);
        if (std::find(visitedNodes.begin(), visitedNodes.end(), allLinks[minIndex].node1) == visitedNodes.end()) {
            visitedNodes.push_back(allLinks[minIndex].node1);
        }
        if (std::find(visitedNodes.begin(), visitedNodes.end(), allLinks[minIndex].node2) == visitedNodes.end()) {
            visitedNodes.push_back(allLinks[minIndex].node2);
        }
        allLinks.erase(allLinks.begin() + minIndex);
    }

    // 根据最小生成树的链路来安装点到点链路设备
    ns3::PointToPointHelper pointToPoint;
    pointToPoint.SetDeviceAttribute("DataRate", ns3::StringValue("5Mbps"));
    pointToPoint.SetChannelAttribute("Delay", ns3::StringValue("2ms"));
    ns3::NetDeviceContainer devices;
    for (const auto& link : mstLinks) {
        ns3::NodeContainer linkNodes;
        linkNodes.Add(link.node1);
        linkNodes.Add(link.node2);
        ns3::NetDeviceContainer linkDevices = pointToPoint.Install(linkNodes);
        devices.Add(linkDevices);
    }

    // 给所有节点安装网络层协议栈
    ns3::InternetStackHelper stack;
    stack.Install(nodes);

    // 配置IP地址分配
    ns3::Ipv4AddressHelper address;
    address.SetBase("10.1.1.0", "255.255.255.0");
    ns3::Ipv4InterfaceContainer interfaces = address.Assign(devices);

    // 在某个节点（这里选择节点0）上安装UDP回声服务器应用
    ns3::UdpEchoServerHelper echoServer(9);
    ns3::ApplicationContainer serverApps = echoServer.Install(nodes.Get(0));
    serverApps.Start(ns3::Seconds(1.0));
    serverApps.Stop(ns3::Seconds(10.0));

    // 在另一个节点（比如节点1）上安装UDP回声客户端应用，指向服务器所在节点的IP和端口
    ns3::UdpEchoClientHelper echoClient(interfaces.GetAddress(0), 9);
    echoClient.SetAttribute("MaxPackets", ns3::UintegerValue(1));
    echoClient.SetAttribute("Interval", ns3::TimeValue(ns3::Seconds(1.0)));
    echoClient.SetAttribute("PacketSize", ns3::UintegerValue(1024));
    ns3::ApplicationContainer clientApps = echoClient.Install(nodes.Get(1));
    clientApps.Start(ns3::Seconds(2.0));
    clientApps.Stop(ns3::Seconds(9.0));

    // 运行仿真
    ns3::Simulator::Run();
    ns3::Simulator::Destroy();

    return 0;
}