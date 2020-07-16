# -*- coding: utf-8 -*-

'''
@Time    : 2020/7/15 14:57
@Author  : panda
@Email   ：panda_isat@yeah.net
@FileName: H_scale.py.py
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
    brdf_file='H:\\BRDF\\dunhuang\\2019年试验\\0824-16.43定点不同高度-场外\\16.43定点不同高度-场外.csv'
    df = pd.read_csv(brdf_file, encoding='utf-8-sig')
    wave_num=800
    point_wave=df.loc[wave_num+5].astype("float")
    wave= pd.array(df['高度'][5:].astype("float"))
    time_ymd=df.loc[1][2:-1]
    time_hms = df.loc[2][2:-1]
    #print(time_ymd,time_hms)
    lat =df.loc[3][2:-1].astype("float")# 纬度
    lon =df.loc[4][2:-1].astype("float")# 经度
    #print('UAV_wavelength:',wave[wave_num])
    ref=np.array(point_wave[2:-1])
    print('ref:', len(ref))
    #height = [10, 10, 20, 20, 30, 30, 50, 50, 100, 100, 150, 150, 200, 200, 250, 250, 250, 250, 200, 200, 150, 150, 100,
    #          100, 50, 50, 30, 30, 20, 20, 10, 10]
    height=['10','20','30','50','100','150','200','250','降250','降200','降150','降100','降50','降30','降20','降10']
    ref=np.array(df.loc[wave_num+5].astype("float")[2:-1])
    print(wave[wave_num],ref)
    h_ref=[]
    for i in range(round(len(ref)/2)):
        #time_ymd = df.loc[0][i*2]
        mean_ref=(ref[i*2]+ref[i*2+1])/2
        h_ref.append(mean_ref)
        print(i,height[i],mean_ref)
    fig, ax = plt.subplots()
    #print(len(local_point_num), len(ref))
    #line1, = ax.plot(local_point_num, local_ref, 'o', ls='-',label=u'跑场，  其归一化std为:'+str(local_std))
    #line1.set_dashes([2, 2, 10, 2])  # 2pt line, 2pt break, 10pt line, 2pt break
    line2, = ax.plot(height, h_ref, 'o', ls='-', dashes=[6, 2], label=u'外场')
    ax.set_xlabel("离地高度（m）")
    ax.set_title("场地均一性，波长："+str(wave[wave_num])+'nm')
    ax.set_ylabel("反射率（%）")

    ax.legend()
    plt.show()
    spice.kclear()
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