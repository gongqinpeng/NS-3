#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/mobility-module.h"
#include "ns3/wifi-module.h"
#include "ns3/aodv-module.h"
#include "ns3/applications-module.h"
#include "ns3/flow-monitor-module.h"

#include <iostream>
#include <fstream>

using namespace ns3;

NS_LOG_COMPONENT_DEFINE("AodvMobilityExample");

void PrintNodePositions(NodeContainer nodes)
{
    for (uint32_t i = 0; i < nodes.GetN(); ++i)
    {
        Ptr<MobilityModel> mobility = nodes.Get(i)->GetObject<MobilityModel>();
        Vector pos = mobility->GetPosition();
        std::cout << "Node " << i << " at (" << pos.x << ", " << pos.y << ")" << std::endl;
    }
}

int main(int argc, char *argv[])
{
    // 定义一些参数
    uint32_t nNodes = 10; // 节点数量
    double simulationTime = 400.0; // 仿真时间，单位为秒
    double dataRate = 256000; // 应用数据速率，单位为bps
    uint32_t packetSize = 512; // 数据包大小，单位为字节
    double interval = 0.1; // 数据包发送间隔，单位为秒
    double maxRange = 70.0; // 通信范围，单位为米

    // 解析命令行参数
    CommandLine cmd(__FILE__);
    cmd.AddValue("nNodes", "Number of nodes", nNodes);
    cmd.AddValue("time", "Simulation time, s", simulationTime);
    cmd.AddValue("maxRange", "Maximum communication range, m", maxRange);
    cmd.Parse(argc, argv);

    // 创建节点
    NodeContainer nodes;
    nodes.Create(nNodes);

    // 配置Wi-Fi
    WifiHelper wifi;
    wifi.SetStandard(WIFI_STANDARD_80211g); // 改为802.11g，提供更高的速率
    YansWifiPhyHelper wifiPhy;
    YansWifiChannelHelper wifiChannel;
    wifiChannel.SetPropagationDelay("ns3::ConstantSpeedPropagationDelayModel");
    wifiChannel.AddPropagationLoss("ns3::RangePropagationLossModel", "MaxRange", DoubleValue(maxRange)); // 设置最大通信范围
    wifiPhy.SetChannel(wifiChannel.Create());
    WifiMacHelper wifiMac;
    wifiMac.SetType("ns3::AdhocWifiMac");
    NetDeviceContainer devices = wifi.Install(wifiPhy, wifiMac, nodes);

    // 配置AODV路由
    AodvHelper aodv;
    InternetStackHelper internet;
    internet.SetRoutingHelper(aodv);
    internet.Install(nodes);

    // 配置IP地址
    Ipv4AddressHelper ipv4;
    Ipv4InterfaceContainer interfaces;
    for (uint32_t i = 0; i < nNodes; ++i)
    {
        std::ostringstream subnet;
        subnet << "10.1." << i + 1 << ".0";
        ipv4.SetBase(subnet.str().c_str(), "255.255.255.0");
        Ipv4InterfaceContainer iface = ipv4.Assign(devices.Get(i));
        interfaces.Add(iface);
    }

    // 配置节点移动模型
    MobilityHelper mobility;
    mobility.SetPositionAllocator("ns3::RandomRectanglePositionAllocator",
                                  "X", StringValue("ns3::UniformRandomVariable[Min=0.0|Max=100.0]"),
                                  "Y", StringValue("ns3::UniformRandomVariable[Min=0.0|Max=100.0]"));
    mobility.SetMobilityModel("ns3::RandomWalk2dMobilityModel",
                              "Mode", StringValue("Time"),
                              "Time", StringValue("20s"), // 节点每20秒移动一次
                              "Speed", StringValue("ns3::UniformRandomVariable[Min=0.5|Max=1.0]"),
                              "Bounds", RectangleValue(Rectangle(0, 100, 0, 100))); // 减少移动速度和范围
    mobility.Install(nodes);

    // 配置流量生成器
    for (uint32_t i = 0; i < nNodes; ++i)
    {
        for (uint32_t j = i + 1; j < nNodes; ++j)
        {
            Ptr<MobilityModel> mobilityI = nodes.Get(i)->GetObject<MobilityModel>();
            Ptr<MobilityModel> mobilityJ = nodes.Get(j)->GetObject<MobilityModel>();

            if (mobilityI->GetDistanceFrom(mobilityJ) <= maxRange) // 检查是否在通信范围内
            {
                uint16_t port = 9000 + i; // 确保每个节点绑定到不同的端口
                OnOffHelper onoff("ns3::UdpSocketFactory", Address(InetSocketAddress(interfaces.GetAddress(j), port)));
                onoff.SetConstantRate(DataRate(dataRate), packetSize);
                ApplicationContainer apps = onoff.Install(nodes.Get(i));
                apps.Start(Seconds(1.0));
                apps.Stop(Seconds(simulationTime - 1.0));

                PacketSinkHelper sink("ns3::UdpSocketFactory", InetSocketAddress(interfaces.GetAddress(i), port));
                apps = sink.Install(nodes.Get(j));
                apps.Start(Seconds(0.0));
            }
        }
    }

    // 启用详细日志记录
    LogComponentEnable("AodvRoutingProtocol", LOG_LEVEL_INFO);
    LogComponentEnable("UdpSocketImpl", LOG_LEVEL_INFO);
    LogComponentEnable("OnOffApplication", LOG_LEVEL_INFO); // 启用流量生成器的日志

    // 启用Flow Monitor
    FlowMonitorHelper flowmon;
    Ptr<FlowMonitor> monitor = flowmon.InstallAll();

    // 在仿真开始时和结束时调用PrintNodePositions
    Simulator::Schedule(Seconds(1.0), &PrintNodePositions, nodes);
    Simulator::Schedule(Seconds(simulationTime - 1.0), &PrintNodePositions, nodes);

    Simulator::Stop(Seconds(simulationTime));
    Simulator::Run();

    // 收集和输出统计数据
    monitor->CheckForLostPackets();
    Ptr<Ipv4FlowClassifier> classifier = DynamicCast<Ipv4FlowClassifier>(flowmon.GetClassifier());
    std::map<FlowId, FlowMonitor::FlowStats> stats = monitor->GetFlowStats();

    for (std::map<FlowId, FlowMonitor::FlowStats>::const_iterator i = stats.begin(); i != stats.end(); ++i)
    {
        Ipv4FlowClassifier::FiveTuple t = classifier->FindFlow(i->first);
        std::cout << "Flow " << i->first << " (" << t.sourceAddress << " -> " << t.destinationAddress << ")\n";
        std::cout << "  Tx Bytes: " << i->second.txBytes << "\n";
        std::cout << "  Rx Bytes: " << i->second.rxBytes << "\n";

        // 改进吞吐量计算
        double duration = i->second.timeLastRxPacket.GetSeconds() - i->second.timeFirstTxPacket.GetSeconds();
        if (duration > 0)
        {
            double throughput = i->second.rxBytes * 8.0 / duration / 1000 / 1000; // Mbps
            std::cout << "  Throughput: " << throughput << " Mbps\n";
        }
        else
        {
            std::cout << "  Throughput: N/A\n";
        }

        // 增加时延统计信息
        if (i->second.rxPackets > 0)
        {
            double avgDelay = i->second.delaySum.GetSeconds() / i->second.rxPackets;
            std::cout << "  Average Delay: " << avgDelay << " s\n";
            std::cout << "  Lost Packets: " << i->second.lostPackets << "\n";
        }
        else
        {
            std::cout << "  Delay: N/A\n";
        }
    }

    Simulator::Destroy();
    return 0;
}
