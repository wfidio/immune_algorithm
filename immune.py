# -*- coding: UTF-8 -*-

######免疫算法求解 TSP 问题##########
import math
import random ###python 绘图模块
import pylab as pl
#参数配置
numOfCity = 50 ##城市数目
maxValue = 100 ##权值最大值
valList = [] ##储存各城市之间路径权值的数组
sizeOfAntiBodyList = 20 ##抗体群大小
iterateTime = 1000 ##迭代次数
avgAffinityList = []
maxAffinityList = []


# ##抗体促进与抑制步骤扩增公式的参数α,β,γ
paraA = 0.9
paraB = 0.5
paraC = 0.4

##随机生成权值图
for row in range(0,numOfCity):
    valList.append([])
    for col in range(0,row+1):
        if(row==col):
            valList[row].append(0)
        else:
            valList[row].append(random.randint(1,maxValue))


##取编号为 x，y 城市之间路径的权值
def val(x,y):
    if(x>y):
        return valList[x-1][y-1]
    elif(x==y):
        return 0
    elif(x<y):
        return val(y,x)


##打印出权值图
graphStr = ''
for i in range(1, numOfCity):
    for j in range(1, numOfCity):
        graphStr = graphStr + str(val(i, j)) + '\t'
    graphStr = graphStr + '\n'
print(graphStr)


##随机生成抗体群,抗体编码从 1 开始
def produceAntibodylist():
    antibodyList = []
    for n in range(0,sizeOfAntiBodyList):
        tempList = []
        for i in range(0, numOfCity):
            tempList.append(i + 1)
        resList = [1]
        tempList.remove(1)
        count = len(tempList)
        for i in range(0,numOfCity-1):
            pos = random.randint(0,count-1)
            resList.append(tempList[pos])
            tempList.remove(tempList[pos])
            count = count-1
        antibodyList.append(resList)
    return antibodyList

##亲和力计算
def affinity (antibody):
    res = 0
    for i in range(0,len(antibody)-1):
        res = res + val(antibody[i],antibody[i+1])

    res = res + val(antibody[len(antibody)-1],antibody[0])
    res = float(res)
    return round(len(antibody)/res,3)

##抗体间相似度计算
def simularity(antibody1,antibody2):
    threshold = 15.0 ##相似度计算中的距离阈值
    distance = 0 ##抗体之间的距离
    for i in range(0,numOfCity):
        distance = distance + math.pow(antibody1[i]-antibody2[i],2)
        distance = math.sqrt(distance)
    # print distance
    if(distance<=threshold):
        return 1
    else:
        return 0


##浓度计算
def density(antibody,antibodyList):
    res = 0
    for antibody2 in antibodyList:
        res = res + simularity(antibody,antibody2)

    res = float(res)/sizeOfAntiBodyList
    # print res
    return res

##根据抗体群产生抗体信息群
def produceAntibodyInfoList(antibodyList):
    antibodyInfoList = []
    for antibody in antibodyList:
        info = {}
        info['affinity'] = affinity(antibody)
        info['antibody'] = antibody
        info['density'] = density(antibody,antibodyList)
        antibodyInfoList.append(info)
    ##根据亲和力对抗体信息群进行降序排列
    antibodyInfoList.sort(reverse=1,key=lambda antibody:antibody['affinity'])
    return antibodyInfoList

##字符互换免疫算子
def patternChange(antibody):
    i=0
    j=0
    while(i==j):
        i = random.randint(1,numOfCity-1)
        j = random.randint(1,numOfCity-1)

    tempAnitbody = antibody
    temp = tempAnitbody[i]
    tempAnitbody[i] = tempAnitbody[j]
    tempAnitbody[j] = temp
    return tempAnitbody


##抗体促进与抑制机制中计算抗体增长个数
def numOfAntibodyInc(antibody,antibodyList):
    res = math.exp(-1*paraA*affinity(antibody))/(paraB*density(antibody,antibodyList)+paraC)
    res = int(res)
    return res


##更新抗体群
def updateAntibodyList(antibodyList,antibodyInfoList):
    increaseAntibodyList = []
    increaseAntibodyInfoList = []

    for antibodyInfo in antibodyInfoList:
        increaseNum = numOfAntibodyInc(antibodyInfo['antibody'],antibodyList)
        ##抗体增长过程
        for i in range(0,increaseNum):
            increaseAntibodyList.append(antibodyInfo['antibody'])
        ##抗体信息表同时更新
        info = {}
        info['antibody'] = antibodyInfo['antibody']
        info['affinity'] = antibodyInfo['affinity']
        increaseAntibodyInfoList.append(info)
    for antibody in increaseAntibodyList:
        antibodyList.append(antibody)
    for antibodyInfo in increaseAntibodyInfoList:
        antibodyInfoList.append(antibodyInfo)



##程序入口
if __name__ =='__main__':

    ##产生权值矩阵
    # produceValMatrix()
    ##生成初始抗体群
    antibodyList = produceAntibodylist()
    ##迭代计数
    count = 0
    ##根据抗体群生成抗体信息表，并按亲和力大小降序排列
    antibodyInfoList = produceAntibodyInfoList(antibodyList)

    ##开始迭代
    for i in range(0,iterateTime):

        ##抗体群促进抑制过程
        updateAntibodyList(antibodyList,antibodyInfoList)
        ##变异算子执行过程
        mutateAntibodyInfoList = []


        ##计算亲和力平均值
        avgAffinity = 0
        for antibodyInfo in antibodyInfoList:
            avgAffinity = avgAffinity + antibodyInfo['affinity']
        avgAffinity = avgAffinity/len(antibodyInfoList)
        avgAffinityList.append(avgAffinity)


        ##遍历抗体群，按照一定概率进行变异
        for antibodyInfo in antibodyInfoList:
            ##根据亲和力得出变异概率,变异概率的基值为 0.02
            p_mutate = math.pow(0.185 - (antibodyInfo['affinity'] - avgAffinity) / avgAffinity, 2)
            # print p_mutate
            if (random.randrange(0, 1000) / 1000.0 < p_mutate):
                mutateAntibody = patternChange(antibodyInfo['antibody'])
                antibodyInfoList.remove(antibodyInfo)
                info = {}
                info['antibody'] = mutateAntibody
                info['affinity'] = affinity(mutateAntibody)
                mutateAntibodyInfoList.append(info)

        for antibodyInfo in mutateAntibodyInfoList:
            antibodyInfoList.append(antibodyInfo)
            ##对变异后的抗体群按照亲和力大小降序排列
        antibodyInfoList.sort(reverse=1, key=lambda antibody: antibody['affinity'])


        ##形成新的抗体群
        antibodyList = []
        count = 0
        count_2 = 0
        while(count<sizeOfAntiBodyList):
            ##剔除重复抗体
            if(antibodyInfoList[count_2]['antibody'] not in antibodyList):
                antibodyList.append(antibodyInfoList[count_2]['antibody'])
                count = count + 1
            count_2 = count_2 + 1

        antibodyInfoList = produceAntibodyInfoList(antibodyList)
        print antibodyInfoList[0]['affinity']
        print antibodyInfoList

    ###对获取结果进行绘制
    timer = []
    for i in range(0,iterateTime):
        timer.append(i+1)
    pl.figure(1)
    pl.plot(timer,avgAffinityList)
    pl.title('Average Affinity')
        # pl.figure(2)
        # pl.plot(timer,maxAffinityList)
        # pl.title('Max Affinity')
    pl.show()
