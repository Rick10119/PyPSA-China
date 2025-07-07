
import pandas as pd
import numpy as np

# 读取生产参数数据 
df_pp = pd.read_csv("data/loadshedding/production_parameters_forLS_modified.csv")
# 读入产品关系矩阵
df_pm = pd.read_csv("data/loadshedding/production_relationship_matrix_modified.csv")

# 获取生产需求量
demand = df_pp['demand'].tolist()
# 获取产品关系矩阵
production_relationship_matrix = df_pm.iloc[:, 1:].values.tolist()

# 向量对应元素乘以矩阵对应列
production_demand_matrix = np.multiply(demand, production_relationship_matrix)

# production_demand_matrix每一行求和(跳过对角元)，计算每个产品的生产需求量
production_demand = np.sum(production_demand_matrix, axis=1) - np.diag(production_demand_matrix)

# 计算每个产品的净需求
net_demand = demand - production_demand
# 计算每个产品的净需求率
net_demand_rate = np.where(np.array(demand) == 0, 0, net_demand / np.array(demand))

# 打印净需求率结果
print("Net Demand Rate:")
for i, rate in enumerate(net_demand_rate):
    print(f"Product {df_pp['materials'].iloc[i]}: {rate:.2%}")

# 结果写入csv文件的对应列
df_pp['net_demand_rate'] = net_demand_rate
df_pp.to_csv("data/loadshedding/production_parameters_forLS_modified.csv", index=False)

# 写完后打印确认信息
print("Net demand rate calculation completed and results saved to CSV.")