#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/mobility-module.h"
#include "ns3/wifi-module.h"
#include "ns3/olsr-helper.h"
#include "ns3/applications-module.h"
#include "ns3/flow-monitor-module.h"

#include <iostream>
#include <fstream>

using namespace ns3;

NS_LOG_COMPONENT_DEFINE("OlsrMobilityExample");

int main(int argc, char *argv[])
{
    // 定义一些参数
    uint32_t nNodes = 5;
    double simulationTime = 50.0; // 仿真时间，单位为秒
    double dataRate = 512000; // 应用数据速率，单位为bps
    uint32_t packetSize = 1024; // 数据包大小，单位为字节
    double interval = 0.1; // 数据包发送间隔，单位为秒

    // 解析命令行参数
    CommandLine cmd(__FILE__);
    cmd.AddValue("nNodes", "Number of nodes", nNodes);
    cmd.AddValue("time", "Simulation time, s", simulationTime);
    cmd.Parse(argc, argv);

    // 创建节点
    NodeContainer nodes;
    nodes.Create(nNodes);

    // 配置Wi-Fi
    WifiHelper wifi;
    wifi.SetStandard(WIFI_STANDARD_80211b);
    YansWifiPhyHelper wifiPhy;
    YansWifiChannelHelper wifiChannel;
    wifiChannel.SetPropagationDelay("ns3::ConstantSpeedPropagationDelayModel");
    wifiChannel.AddPropagationLoss("ns3::TwoRayGroundPropagationLossModel");
    wifiPhy.SetChannel(wifiChannel.Create());
    wifiPhy.Set("TxPowerStart", DoubleValue(16.0206));
    wifiPhy.Set("TxPowerEnd", DoubleValue(16.0206));
    WifiMacHelper wifiMac;
    wifiMac.SetType("ns3::AdhocWifiMac");
    NetDeviceContainer devices = wifi.Install(wifiPhy, wifiMac, nodes);

    // 配置OLSR路由
    OlsrHelper olsr;
    Ipv4ListRoutingHelper list;
    Ipv4StaticRoutingHelper staticRouting;

    list.Add(staticRouting, 0);
    list.Add(olsr, 10);

    InternetStackHelper internet;
    internet.SetRoutingHelper(list);
    internet.Install(nodes);

    // 配置IP地址
    Ipv4AddressHelper ipv4;
    ipv4.SetBase("10.1.1.0", "255.255.255.0");
    Ipv4InterfaceContainer interfaces = ipv4.Assign(devices);

    // 配置节点移动模型
    MobilityHelper mobility;
    mobility.SetPositionAllocator("ns3::RandomRectanglePositionAllocator",
                                  "X", StringValue("ns3::UniformRandomVariable[Min=0.0|Max=500.0]"),
                                  "Y", StringValue("ns3::UniformRandomVariable[Min=0.0|Max=500.0]"));
    mobility.SetMobilityModel("ns3::RandomWalk2dMobilityModel",
                              "Mode", StringValue("Time"),
                              "Time", StringValue("2s"),
                              "Speed", StringValue("ns3::ConstantRandomVariable[Constant=1.0]"),
                              "Bounds", RectangleValue(Rectangle(-50, 550, -50, 550)));
    mobility.Install(nodes);

    // 配置流量生成器
    OnOffHelper onoff("ns3::UdpSocketFactory", Address(InetSocketAddress(interfaces.GetAddress(nNodes-1), 9)));
    onoff.SetConstantRate(DataRate(dataRate), packetSize);
    ApplicationContainer apps = onoff.Install(nodes.Get(0));
    apps.Start(Seconds(0.5));
    apps.Stop(Seconds(simulationTime));

    // 设置接收端
    PacketSinkHelper sink("ns3::UdpSocketFactory", InetSocketAddress(Ipv4Address::GetAny(), 9));
    apps = sink.Install(nodes.Get(nNodes-1));
    apps.Start(Seconds(0.0));

    // 启用Flow Monitor
    FlowMonitorHelper flowmon;
    Ptr<FlowMonitor> monitor = flowmon.InstallAll();

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
        std::cout << "  Throughput: " << i->second.rxBytes * 8.0 / (i->second.timeLastRxPacket.GetSeconds() - i->second.timeFirstTxPacket.GetSeconds()) / 1000 / 1000 << " Mbps\n";
        std::cout << "  Delay: " << i->second.delaySum.GetSeconds() / i->second.rxPackets << " s\n";
    }

    Simulator::Destroy();
    return 0;
}
