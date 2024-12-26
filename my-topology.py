from ns.core import *
from ns.network import *
from ns.internet import *
from ns.point_to_point import *
from ns.applications import *
from ns.internet_apps import *

# 创建节点
nodes = NodeContainer()
nodes.Create(2)  # 创建两个节点

# 设置点对点链路
p2p = PointToPointHelper()
p2p.SetDeviceAttribute("DataRate", StringValue("5Mbps"))
p2p.SetChannelAttribute("Delay", StringValue("2ms"))

# 安装设备
devices = p2p.Install(nodes)

# 设置IP地址
stack = InternetStackHelper()
stack.Install(nodes)
address = Ipv4AddressHelper()
address.SetBase("10.1.1.0", "255.255.255.0")
interfaces = address.Assign(devices)

# 配置应用（例如Ping）
ping = V4PingHelper(interfaces.GetAddress(1))
apps = ping.Install(nodes.Get(0))
apps.Start(Seconds(1.0))
apps.Stop(Seconds(10.0))

# 设置模拟时间和运行
Simulator.Stop(Seconds(11.0))
Simulator.Run()
Simulator.Destroy()

print("Simulation completed")
