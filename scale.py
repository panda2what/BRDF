# -*- coding: utf-8 -*-

'''
@Time    : 2020/7/14 10:27
@Author  : panda
@Email   ：panda_isat@yeah.net
@FileName: scale.py.py
@Software: PyCharm
'''
import numpy as np
import pandas as pd
import spiceypy as spice
import math
import pylab

import ch
ch.set_ch()
import matplotlib.pyplot as plt
def main():
    f_KERNELS = 'D:\\naif\\icy\\kernels\\gen\\'
    spice.furnsh(f_KERNELS + 'spk\\de421.bsp')
    spice.furnsh(f_KERNELS + 'lsk\\naif0011.tls')
    spice.furnsh(f_KERNELS + 'pck\\moon_pa_de421_1900-2050.bpc')
    spice.furnsh(f_KERNELS + 'pck\\earth_latest_high_prec.bpc')
    spice.furnsh(f_KERNELS + 'pck\\pck00010.tpc')
    spice.furnsh(f_KERNELS + 'fk\\moon_080317.tf')
    spice.furnsh(f_KERNELS + 'fk\\moon_assoc_pa.tf')
    spice.furnsh(f_KERNELS + 'fk\\moon_assoc_me.tf')
    spice.furnsh(f_KERNELS + 'lsk\\naif0011.tls.pc')

    #####无人机
    brdf_file='H:\\BRDF\\dunhuang\\2019年试验\\0824-16.13均一性\\16.13均一性.csv'
    #wb_name='16.13均一性'
    #df=pd.read_csv(brdf_file,encoding='utf-8-sig').dropna(axis=1)
    df = pd.read_csv(brdf_file, encoding='utf-8-sig')
    wave_num=200
    point_wave=df.loc[wave_num+5].astype("float")
    wave= pd.array(df['位置'][5:].astype("float"))
    time_ymd=df.loc[1][1:]
    time_hms = df.loc[2][1:]
    lat =df.loc[3][1:].astype("float")# 纬度
    lon =df.loc[4][1:].astype("float")# 经度
    print('UAV_wavelength:',wave[wave_num])
    UAV_ref=np.array(point_wave[2:-1])
    print('UAV_ref:', UAV_ref)
    print('UAV_std:',np.std(UAV_ref,ddof=1))
    UAV_std = round(np.std(UAV_ref, ddof=1), 3)
    for i in range(len(time_ymd)):
        utc=str(time_ymd[i]).replace('/', '-') + 'T' + str(time_hms[i])#2019-8-24T4:25:46未补零
        #et = spice.spiceypy.utc2et(utc)
        solar_azimuth, solar_zenith = C_solar_zenith(utc, lat[i], lon[i])#方位角，高度角
        #print(utc, lat[i], lon[i], solar_azimuth, solar_zenith)#
    #跑场
    brdf_file='H:\\BRDF\\dunhuang\\2019年试验\\0824-17.01与地面做对比\\17.01与地面做对比.csv'
    df = pd.read_csv(brdf_file, encoding='utf-8-sig')

    point_wave=df.loc[wave_num+5].astype("float")
    wave= pd.array(df['位置'][5:].astype("float"))
    time_ymd=df.loc[1][1:]
    time_hms = df.loc[2][1:]
    lat =df.loc[3][1:].astype("float")# 纬度
    lon =df.loc[4][1:].astype("float")# 经度
    print('local_wavelength:',wave[wave_num])
    local_ref=np.array(point_wave[2:-1])
    print('local_ref:', local_ref)
    print('local_std:',np.std(local_ref,ddof=1))
    local_std=round(np.std(local_ref, ddof=1),3)
    for i in range(len(time_ymd)):
        utc=str(time_ymd[i]).replace('/', '-') + 'T' + str(time_hms[i])#2019-8-24T4:25:46未补零
        #et = spice.spiceypy.utc2et(utc)
        solar_azimuth, solar_zenith = C_solar_zenith(utc, lat[i], lon[i])#方位角，高度角
        #print(utc, lat[i], lon[i], solar_azimuth, solar_zenith)#
    spice.kclear()
    BB_file='H:\\BRDF\\dunhuang\\2019年试验\\2019年白板数据\\MFB99-359Y-19.csv'#白板数据
    df = pd.read_csv(BB_file, encoding='utf-8-sig')
    point_wave = df.loc[15].astype("float")
    local_point_num=[1,2,3,4,5,6,7,8,9,10,11,12]
    UAV_point_num = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,13,14,15,16,17,18]
    fig, ax = plt.subplots()
    #print(len(local_point_num), len(ref))
    #line1, = ax.plot(local_point_num, local_ref, 'o', ls='-',label=u'跑场，  其归一化std为:'+str(local_std))
    #line1.set_dashes([2, 2, 10, 2])  # 2pt line, 2pt break, 10pt line, 2pt break
    line2, = ax.plot(UAV_point_num, UAV_ref, 'o', ls='-', dashes=[6, 2], label=u'无人机,其std为:'+str(UAV_std))
    ax.set_xlabel("点编号")
    ax.set_title("场地均一性，波长："+str(wave[wave_num])+'nm')
    ax.set_ylabel("反射率（%）")
    ax.set_xlim(0, 18)
    ax.legend()
    plt.show()

    #print(point_wave)### ax.plot(x, y, 'o', ls='-', ms=4, markevery=case)
    print('##FINISH###')

def C_solar_zenith(utc, lat, lon):
    # H_S太阳高度角，lat:地理纬度，solar_RA太阳赤纬，t时角,lon经度：94
    TimeZone = 8.
    et = spice.spiceypy.utc2et(utc)
    utc = (spice.spiceypy.et2utc(et, 'ISOC', 0, 24))  # 时间转换
    year = float(utc[0:4])
    # print(year)
    et1 = spice.spiceypy.utc2et(utc[0:4] + '-01-01T00:00:00')
    #print( (spice.spiceypy.et2utc(et1, 'ISOC', 0, 24)) ) # 时间转换

    DOY = (math.ceil((et - et1) / 60 / 60 / 24))
    # print("DOY:",DOY)
    N0 = 79.6764 + 0.2422 * (year - 1985) - np.floor((year - 1985) / 4.0)
    # print("N0:",N0)
    sitar = 2 * np.pi * (DOY - N0) / 365.2422
    ED = 0.3723 + 23.2567 * np.sin(sitar) + 0.1149 * np.sin(2 * sitar) - 0.1712 * np.sin(3 * sitar) - 0.758 * np.cos(
        sitar) + 0.3656 * np.cos(2 * sitar) + 0.0201 * np.cos(3 * sitar)
    ED = ED * np.pi / 180.
    if (lon >= 0):
        if (TimeZone == -13):
            dLon = lon - (np.floor((lon * 10. - 75.) / 150.) + 1.) * 15.0
        else:
            dLon = lon - TimeZone * 15.0
    else:
        if (TimeZone == -13):
            dLon = (np.floor((lon * 10. - 75.) / 150.) + 1) * 15.0 - lon
        else:
            dLon = TimeZone * 15.0 - lon
    Et = 0.0028 - 1.9857 * np.sin(sitar) + 9.9059 * np.sin(2 * sitar) - 7.0924 * np.cos(sitar) - 0.6882 * np.cos(
        2 * sitar)

    hour = float(utc[11:13])
    mins = float(utc[14:16])
    sec = float(utc[17:19])

    gtdt = hour + mins / 60.0 + sec / 3600.0 + dLon / 15.#地方时
    # print('gtdt: ',gtdt )
    gtdt = gtdt + Et / 60.0
    dTimeAngle1 = 15.0 * (gtdt-12.)  #########时角和会有区别
    # print('dTimeAngle:',dTimeAngle)
    dTimeAngle = dTimeAngle1 * np.pi / 180.
    latitudeArc = lat * np.pi / 180.
    # print('lat:',latitudeArc,lat)
    HeightAngleArc = np.arcsin(
        np.sin(latitudeArc) * np.sin(ED) + np.cos(latitudeArc) * np.cos(ED) * np.cos(dTimeAngle))  # % //高度角计算公式
    CosAzimuthAngle = (np.sin(HeightAngleArc) * np.sin(latitudeArc) - np.sin(ED)) / np.cos(HeightAngleArc) / np.cos(
        latitudeArc)  # % //方位角计算公式
    AzimuthAngleArc = np.arccos(CosAzimuthAngle)  #
    HeightAngle = 90-HeightAngleArc * 180 / np.pi  #
    AzimuthAngle = AzimuthAngleArc * 180 / np.pi  #
    if (dTimeAngle < 0):
        AzimuthAngle = 180 - AzimuthAngle
    else:
        AzimuthAngle = 180 + AzimuthAngle
    solar_azimuth = AzimuthAngle
    solar_zenith = HeightAngle
    return (solar_azimuth, solar_zenith)



if __name__=='__main__':
    main()
