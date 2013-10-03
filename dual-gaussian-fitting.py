#
#   A program to attempt to fit two gaussians to data, by first fitting one, then fitting a second.
#
#
#                 by Rob Heath
#

#
#   Currently this finds the FIRST PEAK and fits to the TAIL SIDE of this - I really need to find BOTH peaks...
#     ... which means differentiating and looking for sign changes.
#
#
#          ###### This code was abandoned as there are easier routes to success (OriginLab 8.5) #########

import numpy as np
import scipy as sp
import csv

filename = 'output.txt'

y_points = [809,884,854,831,838,798,862,890,829,837,867,839,829,874,868,833,822,869,795,810,802,790,850,839,868,829,871,836,811,850,848,810,822,844,829,861,872,890,872,847,879,817,817,882,782,852,839,874,850,808,910,873,900,850,857,851,889,823,856,871,852,868,903,875,856,902,854,858,888,854,889,856,851,842,909,848,880,864,825,873,898,882,863,851,838,856,899,843,854,805,845,864,897,841,889,844,817,908,841,939,864,884,825,895,868,857,887,847,919,846,873,893,900,884,899,884,877,869,924,856,892,849,891,877,868,858,890,863,917,877,904,902,891,879,900,940,948,880,851,858,900,887,902,903,939,832,852,857,927,897,891,911,897,914,915,831,861,922,887,854,850,926,823,914,916,942,919,908,901,852,908,939,917,884,919,876,961,919,971,975,955,942,929,991,872,957,913,926,901,951,927,936,919,893,915,941,985,983,911,939,1070,938,950,922,948,973,968,951,909,964,982,987,926,1023,971,952,932,936,955,922,974,932,975,960,958,936,957,900,930,906,957,1006,971,943,985,946,957,990,1003,945,943,950,969,1026,962,923,944,987,1003,1000,986,949,1001,977,992,1008,972,1032,1021,1064,990,1005,1040,1051,1024,1057,1042,1068,1044,1052,1128,1058,1125,1084,1114,1086,1139,1189,1166,1185,1212,1140,1215,1188,1158,1238,1208,1287,1285,1338,1212,1339,1357,1390,1424,1379,1424,1506,1534,1514,1568,1581,1635,1689,1632,1659,1796,1820,1780,1808,1905,1891,1973,2031,2086,2167,2085,2285,2311,2348,2399,2486,2537,2637,2648,2659,2815,2916,2994,3033,3238,3261,3296,3449,3487,3574,3578,3744,3816,3842,4011,4108,4274,4264,4461,4327,4579,4651,4696,4656,4852,4970,5026,4970,5310,5111,5413,5283,5462,5359,5723,5651,5848,5813,5808,5711,5894,5962,6130,6042,6157,6152,6109,6031,6301,6285,6335,6276,6301,6343,6339,6390,6510,6336,6609,6452,6533,6490,6404,6433,6448,6573,6450,6566,6393,6522,6569,6591,6724,6608,6540,6680,6541,6593,6757,6797,6783,6781,6850,6849,6912,7078,7150,7179,7223,7370,7637,7764,7994,8154,8317,8529,8623,8842,9197,9473,9638,9814,9823,10184,10459,10546,10956,11178,11447,11536,11688,11741,11860,11944,12078,12107,12194,12332,12430,12307,12420,12325,12181,12142,12285,12082,11961,11890,11803,11661,11678,11578,11252,11323,10949,10708,10563,10461,10204,10166,9821,9749,9581,9434,9379,8888,8828,8681,8536,8449,8090,8005,7912,7857,7656,7323,7259,7032,6869,6833,6747,6481,6510,6314,6166,5999,5859,5806,5729,5639,5381,5387,5276,5119,5088,4984,4870,4831,4815,4532,4616,4436,4255,4278,4151,4151,4019,3999,3865,3739,3766,3744,3702,3531,3431,3455,3377,3399,3323,3147,3121,3266,3151,3066,2915,2985,2823,2834,2849,2843,2772,2602,2645,2588,2646,2597,2587,2543,2535,2452,2367,2368,2353,2193,2231,2347,2259,2143,2241,2219,2174,2151,2201,2060,2145,2089,1986,2000,1988,1986,1923,1984,1868,1956,1888,1865,1980,1863,1886,1925,1861,1819,1801,1862,1785,1791,1713,1825,1790,1769,1788,1730,1745,1754,1708,1820,1697,1676,1707,1698,1759,1787,1718,1677,1706,1651,1643,1753,1715,1638,1657,1738,1653,1679,1665,1693,1596,1660,1713,1669,1638,1649,1614,1673,1617,1623,1575,1614,1547,1669,1645,1682,1638,1607,1625,1588,1577,1684,1622,1583,1621,1569,1593,1642,1562,1589,1588,1619,1621,1592,1616,1599,1533,1620,1590,1591,1500,1590,1655,1585,1589,1540,1575,1655,1557,1585,1649,1499,1531,1647,1520,1575,1623,1590,1584,1564,1524,1602,1579,1514,1566,1525,1466,1572,1603,1559,1515,1624,1532,1526,1452,1541,1527,1500,1473,1548,1480,1507,1518,1539,1487,1482,1423,1437,1444,1406,1503,1407,1518,1407,1499,1522,1452,1381,1455,1512,1391,1363,1357,1454,1436,1367,1480,1420,1427,1395,1416,1354,1442,1333,1378,1396,1295,1332,1301,1304,1337,1350,1353,1348,1311,1299,1262,1343,1317,1251,1249,1304,1258,1299,1220,1360,1192,1229,1287,1247,1218,1243,1262,1211,1247,1257,1221,1199,1217,1216,1182,1175,1182,1251,1137,1214,1158,1158,1198,1215,1201,1151,1208,1112,1130,1187,1184,1133,1145,1191,1126,1058,1120,1093,1083,1106,1071,1039,1069,1057,1108,1070,1028,1051,1165,1122,1008,1086,1132,1040,1082,1064,1100,1043,1008,1034,1099,1035,992,1011,1034,1054,994,998,993,1072,982,1005,970,965,1055,964,999,1020,999,986,949,963,1018,991,953,970,950,974,1057,1008,999,968,927,938,946,978,983,970,954,929,964,982,942,948,946,922,1040,946,981,946,914,959,934,976,939,990,944,952,933,907,923,924,962,939,936,926,930,909,919,960,935,920,890,960,928,896,890,910,876,907,922,928,930,905,907,941,879,926,926,862,850,941,879,913,973,902,873,964,928,897,860,916,951,889,967,994,892,911,940,924,913,909,941,985,932,950,921,960,997,943,951,908,892,960,929,963,1016,971,1037,978,957,940,1025,965,978,911,970,997,974,985,986,1011,1017,997,988,985,1000,980,983,1019,966,1011,1005,1000,965,960,1028,1033,992,995,1059,1017,984,1053,1020,1094,1019,1098,1068,1022,1062,1048,1039,973,1106,1064,1007,1044,1073,1044,998,1061,1020,1076,1096,1200,1119,1179,1103,1146,1128,1163,1086,1126,1145,1117,1110,1092,1119,1158,1166,1160,1153,1143,1127,1142,1221,1126,1170,1155,1179,1171,1214,1125,1228,1194,1186,1200,1210,1233,1246,1189,1262,1163,1236,1216,1220,1260,1277,1193,1271,1178,1264,1227,1228,1202,1247,1293,1253,1282,1240,1297,1280,1241,1250,1264,1325,1302,1251,1239,1275,1298,1273,1317,1332,1226,1295,1326,1323,1305,1339,1291,1322,1351,1396,1342,1352,1328,1392,1373,1342,1343,1322,1407,1328,1413,1367,1339,1342,1354,1317,1372,1366,1429,1379,1371,1433,1468,1464,1379,1404,1392,1428,1445,1440,1424,1444,1423,1436,1397,1500,1473,1472,1445,1509,1453,1503,1474,1479,1471,1528,1459,1543,1489,1524,1498,1514,1512,1483,1549,1502,1548,1472,1493,1524,1465,1586,1508,1490,1528,1517,1541,1503,1541,1564,1509,1577,1451,1601,1584,1545,1484,1551,1557,1541,1494,1594,1483,1566,1556,1561,1658,1601,1518,1587,1560,1588,1624,1494,1562,1580,1509,1572,1543,1528,1613,1580,1505,1576,1509,1522,1595,1540,1555,1569,1482,1581,1518,1583,1576,1536,1528,1524,1520,1560,1535,1574,1532,1512,1497,1522,1524,1546,1544,1566,1500,1581,1601,1562,1541,1553,1504,1573,1543,1499,1554,1529,1541,1549,1489,1605,1530,1568,1526,1577,1573,1505,1584,1542,1529,1514,1541,1478,1490,1487,1524,1561,1515,1467,1473,1483,1546,1498,1536,1552,1542,1505,1520,1491,1451,1493,1456,1439,1432,1432,1423,1460,1445,1501,1453,1479,1537,1474,1447,1508,1476,1406,1427,1466,1427,1432,1459,1451,1475,1491,1388,1448,1383,1429,1417,1468,1374,1459,1376,1428,1406,1344,1287,1411,1354,1401,1394,1421,1348,1411,1371,1353,1290,1393,1296,1346,1388,1402,1344,1378,1381,1339,1318,1319,1399,1389,1353,1352,1324,1358,1384,1294,1315,1310,1267,1320,1276,1269,1328,1330,1308,1284,1231,1315,1272,1269,1238,1278,1299,1290,1289,1261,1277,1273,1243,1223,1228,1314,1209,1250,1271,1197,1218,1150,1162,1185,1182,1166,1284,1204,1178,1221,1221,1129,1126,1195,1197,1213,1200,1135,1109,1150,1128,1169,1149,1136,1155,1095,1116,1120,1125,1061,1141,1088,1079,1107,1135,1074,1073,1066,1021,1104,1075,1043,1009,1016,1075,1019,1030,1118,1051,997,992,1012,1051,1039,1040,1065,1035,1033,1046,1066,1029,1059,988,1010,996,1041,1028,1008,1000,990,1019,1020,920,999,1011,1012,996,976,980,987,959,976,986,951,902,924,974,979,938,893,879,923,951,974,944,898,955,947,898,942,933,903,922,890,868,910,934,823,905,877,916,932,872,890,882,859,888,886,888,906,826,825,904,866,850,853,852,810,853,853,888,820,867,884,774,876,851,823,776,871,870,870,808,846,776,874,814,803,843,806,817,811,789,757]
x_points = [8E-12,1.2E-11,1.6E-11,2E-11,2.4E-11,2.8E-11,3.2E-11,3.6E-11,4E-11,4.4E-11,4.8E-11,5.2E-11,5.6E-11,6E-11,6.4E-11,6.8E-11,7.2E-11,7.6E-11,8E-11,8.4E-11,8.8E-11,9.2E-11,9.6E-11,1E-10,1.04E-10,1.08E-10,1.12E-10,1.16E-10,1.2E-10,1.24E-10,1.28E-10,1.32E-10,1.36E-10,1.4E-10,1.44E-10,1.48E-10,1.52E-10,1.56E-10,1.6E-10,1.64E-10,1.68E-10,1.72E-10,1.76E-10,1.8E-10,1.84E-10,1.88E-10,1.92E-10,1.96E-10,2E-10,2.04E-10,2.08E-10,2.12E-10,2.16E-10,2.2E-10,2.24E-10,2.28E-10,2.32E-10,2.36E-10,2.4E-10,2.44E-10,2.48E-10,2.52E-10,2.56E-10,2.6E-10,2.64E-10,2.68E-10,2.72E-10,2.76E-10,2.8E-10,2.84E-10,2.88E-10,2.92E-10,2.96E-10,3E-10,3.04E-10,3.08E-10,3.12E-10,3.16E-10,3.2E-10,3.24E-10,3.28E-10,3.32E-10,3.36E-10,3.4E-10,3.44E-10,3.48E-10,3.52E-10,3.56E-10,3.6E-10,3.64E-10,3.68E-10,3.72E-10,3.76E-10,3.8E-10,3.84E-10,3.88E-10,3.92E-10,3.96E-10,4E-10,4.04E-10,4.08E-10,4.12E-10,4.16E-10,4.2E-10,4.24E-10,4.28E-10,4.32E-10,4.36E-10,4.4E-10,4.44E-10,4.48E-10,4.52E-10,4.56E-10,4.6E-10,4.64E-10,4.68E-10,4.72E-10,4.76E-10,4.8E-10,4.84E-10,4.88E-10,4.92E-10,4.96E-10,5E-10,5.04E-10,5.08E-10,5.12E-10,5.16E-10,5.2E-10,5.24E-10,5.28E-10,5.32E-10,5.36E-10,5.4E-10,5.44E-10,5.48E-10,5.52E-10,5.56E-10,5.6E-10,5.64E-10,5.68E-10,5.72E-10,5.76E-10,5.8E-10,5.84E-10,5.88E-10,5.92E-10,5.96E-10,6E-10,6.04E-10,6.08E-10,6.12E-10,6.16E-10,6.2E-10,6.24E-10,6.28E-10,6.32E-10,6.36E-10,6.4E-10,6.44E-10,6.48E-10,6.52E-10,6.56E-10,6.6E-10,6.64E-10,6.68E-10,6.72E-10,6.76E-10,6.8E-10,6.84E-10,6.88E-10,6.92E-10,6.96E-10,7E-10,7.04E-10,7.08E-10,7.12E-10,7.16E-10,7.2E-10,7.24E-10,7.28E-10,7.32E-10,7.36E-10,7.4E-10,7.44E-10,7.48E-10,7.52E-10,7.56E-10,7.6E-10,7.64E-10,7.68E-10,7.72E-10,7.76E-10,7.8E-10,7.84E-10,7.88E-10,7.92E-10,7.96E-10,8E-10,8.04E-10,8.08E-10,8.12E-10,8.16E-10,8.2E-10,8.24E-10,8.28E-10,8.32E-10,8.36E-10,8.4E-10,8.44E-10,8.48E-10,8.52E-10,8.56E-10,8.6E-10,8.64E-10,8.68E-10,8.72E-10,8.76E-10,8.8E-10,8.84E-10,8.88E-10,8.92E-10,8.96E-10,9E-10,9.04E-10,9.08E-10,9.12E-10,9.16E-10,9.2E-10,9.24E-10,9.28E-10,9.32E-10,9.36E-10,9.4E-10,9.44E-10,9.48E-10,9.52E-10,9.56E-10,9.6E-10,9.64E-10,9.68E-10,9.72E-10,9.76E-10,9.8E-10,9.84E-10,9.88E-10,9.92E-10,9.96E-10,1E-9,1.004E-9,1.008E-9,1.012E-9,1.016E-9,1.02E-9,1.024E-9,1.028E-9,1.032E-9,1.036E-9,1.04E-9,1.044E-9,1.048E-9,1.052E-9,1.056E-9,1.06E-9,1.064E-9,1.068E-9,1.072E-9,1.076E-9,1.08E-9,1.084E-9,1.088E-9,1.092E-9,1.096E-9,1.1E-9,1.104E-9,1.108E-9,1.112E-9,1.116E-9,1.12E-9,1.124E-9,1.128E-9,1.132E-9,1.136E-9,1.14E-9,1.144E-9,1.148E-9,1.152E-9,1.156E-9,1.16E-9,1.164E-9,1.168E-9,1.172E-9,1.176E-9,1.18E-9,1.184E-9,1.188E-9,1.192E-9,1.196E-9,1.2E-9,1.204E-9,1.208E-9,1.212E-9,1.216E-9,1.22E-9,1.224E-9,1.228E-9,1.232E-9,1.236E-9,1.24E-9,1.244E-9,1.248E-9,1.252E-9,1.256E-9,1.26E-9,1.264E-9,1.268E-9,1.272E-9,1.276E-9,1.28E-9,1.284E-9,1.288E-9,1.292E-9,1.296E-9,1.3E-9,1.304E-9,1.308E-9,1.312E-9,1.316E-9,1.32E-9,1.324E-9,1.328E-9,1.332E-9,1.336E-9,1.34E-9,1.344E-9,1.348E-9,1.352E-9,1.356E-9,1.36E-9,1.364E-9,1.368E-9,1.372E-9,1.376E-9,1.38E-9,1.384E-9,1.388E-9,1.392E-9,1.396E-9,1.4E-9,1.404E-9,1.408E-9,1.412E-9,1.416E-9,1.42E-9,1.424E-9,1.428E-9,1.432E-9,1.436E-9,1.44E-9,1.444E-9,1.448E-9,1.452E-9,1.456E-9,1.46E-9,1.464E-9,1.468E-9,1.472E-9,1.476E-9,1.48E-9,1.484E-9,1.488E-9,1.492E-9,1.496E-9,1.5E-9,1.504E-9,1.508E-9,1.512E-9,1.516E-9,1.52E-9,1.524E-9,1.528E-9,1.532E-9,1.536E-9,1.54E-9,1.544E-9,1.548E-9,1.552E-9,1.556E-9,1.56E-9,1.564E-9,1.568E-9,1.572E-9,1.576E-9,1.58E-9,1.584E-9,1.588E-9,1.592E-9,1.596E-9,1.6E-9,1.604E-9,1.608E-9,1.612E-9,1.616E-9,1.62E-9,1.624E-9,1.628E-9,1.632E-9,1.636E-9,1.64E-9,1.644E-9,1.648E-9,1.652E-9,1.656E-9,1.66E-9,1.664E-9,1.668E-9,1.672E-9,1.676E-9,1.68E-9,1.684E-9,1.688E-9,1.692E-9,1.696E-9,1.7E-9,1.704E-9,1.708E-9,1.712E-9,1.716E-9,1.72E-9,1.724E-9,1.728E-9,1.732E-9,1.736E-9,1.74E-9,1.744E-9,1.748E-9,1.752E-9,1.756E-9,1.76E-9,1.764E-9,1.768E-9,1.772E-9,1.776E-9,1.78E-9,1.784E-9,1.788E-9,1.792E-9,1.796E-9,1.8E-9,1.804E-9,1.808E-9,1.812E-9,1.816E-9,1.82E-9,1.824E-9,1.828E-9,1.832E-9,1.836E-9,1.84E-9,1.844E-9,1.848E-9,1.852E-9,1.856E-9,1.86E-9,1.864E-9,1.868E-9,1.872E-9,1.876E-9,1.88E-9,1.884E-9,1.888E-9,1.892E-9,1.896E-9,1.9E-9,1.904E-9,1.908E-9,1.912E-9,1.916E-9,1.92E-9,1.924E-9,1.928E-9,1.932E-9,1.936E-9,1.94E-9,1.944E-9,1.948E-9,1.952E-9,1.956E-9,1.96E-9,1.964E-9,1.968E-9,1.972E-9,1.976E-9,1.98E-9,1.984E-9,1.988E-9,1.992E-9,1.996E-9,2E-9,2.004E-9,2.008E-9,2.012E-9,2.016E-9,2.02E-9,2.024E-9,2.028E-9,2.032E-9,2.036E-9,2.04E-9,2.044E-9,2.048E-9,2.052E-9,2.056E-9,2.06E-9,2.064E-9,2.068E-9,2.072E-9,2.076E-9,2.08E-9,2.084E-9,2.088E-9,2.092E-9,2.096E-9,2.1E-9,2.104E-9,2.108E-9,2.112E-9,2.116E-9,2.12E-9,2.124E-9,2.128E-9,2.132E-9,2.136E-9,2.14E-9,2.144E-9,2.148E-9,2.152E-9,2.156E-9,2.16E-9,2.164E-9,2.168E-9,2.172E-9,2.176E-9,2.18E-9,2.184E-9,2.188E-9,2.192E-9,2.196E-9,2.2E-9,2.204E-9,2.208E-9,2.212E-9,2.216E-9,2.22E-9,2.224E-9,2.228E-9,2.232E-9,2.236E-9,2.24E-9,2.244E-9,2.248E-9,2.252E-9,2.256E-9,2.26E-9,2.264E-9,2.268E-9,2.272E-9,2.276E-9,2.28E-9,2.284E-9,2.288E-9,2.292E-9,2.296E-9,2.3E-9,2.304E-9,2.308E-9,2.312E-9,2.316E-9,2.32E-9,2.324E-9,2.328E-9,2.332E-9,2.336E-9,2.34E-9,2.344E-9,2.348E-9,2.352E-9,2.356E-9,2.36E-9,2.364E-9,2.368E-9,2.372E-9,2.376E-9,2.38E-9,2.384E-9,2.388E-9,2.392E-9,2.396E-9,2.4E-9,2.404E-9,2.408E-9,2.412E-9,2.416E-9,2.42E-9,2.424E-9,2.428E-9,2.432E-9,2.436E-9,2.44E-9,2.444E-9,2.448E-9,2.452E-9,2.456E-9,2.46E-9,2.464E-9,2.468E-9,2.472E-9,2.476E-9,2.48E-9,2.484E-9,2.488E-9,2.492E-9,2.496E-9,2.5E-9,2.504E-9,2.508E-9,2.512E-9,2.516E-9,2.52E-9,2.524E-9,2.528E-9,2.532E-9,2.536E-9,2.54E-9,2.544E-9,2.548E-9,2.552E-9,2.556E-9,2.56E-9,2.564E-9,2.568E-9,2.572E-9,2.576E-9,2.58E-9,2.584E-9,2.588E-9,2.592E-9,2.596E-9,2.6E-9,2.604E-9,2.608E-9,2.612E-9,2.616E-9,2.62E-9,2.624E-9,2.628E-9,2.632E-9,2.636E-9,2.64E-9,2.644E-9,2.648E-9,2.652E-9,2.656E-9,2.66E-9,2.664E-9,2.668E-9,2.672E-9,2.676E-9,2.68E-9,2.684E-9,2.688E-9,2.692E-9,2.696E-9,2.7E-9,2.704E-9,2.708E-9,2.712E-9,2.716E-9,2.72E-9,2.724E-9,2.728E-9,2.732E-9,2.736E-9,2.74E-9,2.744E-9,2.748E-9,2.752E-9,2.756E-9,2.76E-9,2.764E-9,2.768E-9,2.772E-9,2.776E-9,2.78E-9,2.784E-9,2.788E-9,2.792E-9,2.796E-9,2.8E-9,2.804E-9,2.808E-9,2.812E-9,2.816E-9,2.82E-9,2.824E-9,2.828E-9,2.832E-9,2.836E-9,2.84E-9,2.844E-9,2.848E-9,2.852E-9,2.856E-9,2.86E-9,2.864E-9,2.868E-9,2.872E-9,2.876E-9,2.88E-9,2.884E-9,2.888E-9,2.892E-9,2.896E-9,2.9E-9,2.904E-9,2.908E-9,2.912E-9,2.916E-9,2.92E-9,2.924E-9,2.928E-9,2.932E-9,2.936E-9,2.94E-9,2.944E-9,2.948E-9,2.952E-9,2.956E-9,2.96E-9,2.964E-9,2.968E-9,2.972E-9,2.976E-9,2.98E-9,2.984E-9,2.988E-9,2.992E-9,2.996E-9,3E-9,3.004E-9,3.008E-9,3.012E-9,3.016E-9,3.02E-9,3.024E-9,3.028E-9,3.032E-9,3.036E-9,3.04E-9,3.044E-9,3.048E-9,3.052E-9,3.056E-9,3.06E-9,3.064E-9,3.068E-9,3.072E-9,3.076E-9,3.08E-9,3.084E-9,3.088E-9,3.092E-9,3.096E-9,3.1E-9,3.104E-9,3.108E-9,3.112E-9,3.116E-9,3.12E-9,3.124E-9,3.128E-9,3.132E-9,3.136E-9,3.14E-9,3.144E-9,3.148E-9,3.152E-9,3.156E-9,3.16E-9,3.164E-9,3.168E-9,3.172E-9,3.176E-9,3.18E-9,3.184E-9,3.188E-9,3.192E-9,3.196E-9,3.2E-9,3.204E-9,3.208E-9,3.212E-9,3.216E-9,3.22E-9,3.224E-9,3.228E-9,3.232E-9,3.236E-9,3.24E-9,3.244E-9,3.248E-9,3.252E-9,3.256E-9,3.26E-9,3.264E-9,3.268E-9,3.272E-9,3.276E-9,3.28E-9,3.284E-9,3.288E-9,3.292E-9,3.296E-9,3.3E-9,3.304E-9,3.308E-9,3.312E-9,3.316E-9,3.32E-9,3.324E-9,3.328E-9,3.332E-9,3.336E-9,3.34E-9,3.344E-9,3.348E-9,3.352E-9,3.356E-9,3.36E-9,3.364E-9,3.368E-9,3.372E-9,3.376E-9,3.38E-9,3.384E-9,3.388E-9,3.392E-9,3.396E-9,3.4E-9,3.404E-9,3.408E-9,3.412E-9,3.416E-9,3.42E-9,3.424E-9,3.428E-9,3.432E-9,3.436E-9,3.44E-9,3.444E-9,3.448E-9,3.452E-9,3.456E-9,3.46E-9,3.464E-9,3.468E-9,3.472E-9,3.476E-9,3.48E-9,3.484E-9,3.488E-9,3.492E-9,3.496E-9,3.5E-9,3.504E-9,3.508E-9,3.512E-9,3.516E-9,3.52E-9,3.524E-9,3.528E-9,3.532E-9,3.536E-9,3.54E-9,3.544E-9,3.548E-9,3.552E-9,3.556E-9,3.56E-9,3.564E-9,3.568E-9,3.572E-9,3.576E-9,3.58E-9,3.584E-9,3.588E-9,3.592E-9,3.596E-9,3.6E-9,3.604E-9,3.608E-9,3.612E-9,3.616E-9,3.62E-9,3.624E-9,3.628E-9,3.632E-9,3.636E-9,3.64E-9,3.644E-9,3.648E-9,3.652E-9,3.656E-9,3.66E-9,3.664E-9,3.668E-9,3.672E-9,3.676E-9,3.68E-9,3.684E-9,3.688E-9,3.692E-9,3.696E-9,3.7E-9,3.704E-9,3.708E-9,3.712E-9,3.716E-9,3.72E-9,3.724E-9,3.728E-9,3.732E-9,3.736E-9,3.74E-9,3.744E-9,3.748E-9,3.752E-9,3.756E-9,3.76E-9,3.764E-9,3.768E-9,3.772E-9,3.776E-9,3.78E-9,3.784E-9,3.788E-9,3.792E-9,3.796E-9,3.8E-9,3.804E-9,3.808E-9,3.812E-9,3.816E-9,3.82E-9,3.824E-9,3.828E-9,3.832E-9,3.836E-9,3.84E-9,3.844E-9,3.848E-9,3.852E-9,3.856E-9,3.86E-9,3.864E-9,3.868E-9,3.872E-9,3.876E-9,3.88E-9,3.884E-9,3.888E-9,3.892E-9,3.896E-9,3.9E-9,3.904E-9,3.908E-9,3.912E-9,3.916E-9,3.92E-9,3.924E-9,3.928E-9,3.932E-9,3.936E-9,3.94E-9,3.944E-9,3.948E-9,3.952E-9,3.956E-9,3.96E-9,3.964E-9,3.968E-9,3.972E-9,3.976E-9,3.98E-9,3.984E-9,3.988E-9,3.992E-9,3.996E-9,4E-9,4.004E-9,4.008E-9,4.012E-9,4.016E-9,4.02E-9,4.024E-9,4.028E-9,4.032E-9,4.036E-9,4.04E-9,4.044E-9,4.048E-9,4.052E-9,4.056E-9,4.06E-9,4.064E-9,4.068E-9,4.072E-9,4.076E-9,4.08E-9,4.084E-9,4.088E-9,4.092E-9,4.096E-9,4.1E-9,4.104E-9,4.108E-9,4.112E-9,4.116E-9,4.12E-9,4.124E-9,4.128E-9,4.132E-9,4.136E-9,4.14E-9,4.144E-9,4.148E-9,4.152E-9,4.156E-9,4.16E-9,4.164E-9,4.168E-9,4.172E-9,4.176E-9,4.18E-9,4.184E-9,4.188E-9,4.192E-9,4.196E-9,4.2E-9,4.204E-9,4.208E-9,4.212E-9,4.216E-9,4.22E-9,4.224E-9,4.228E-9,4.232E-9,4.236E-9,4.24E-9,4.244E-9,4.248E-9,4.252E-9,4.256E-9,4.26E-9,4.264E-9,4.268E-9,4.272E-9,4.276E-9,4.28E-9,4.284E-9,4.288E-9,4.292E-9,4.296E-9,4.3E-9,4.304E-9,4.308E-9,4.312E-9,4.316E-9,4.32E-9,4.324E-9,4.328E-9,4.332E-9,4.336E-9,4.34E-9,4.344E-9,4.348E-9,4.352E-9,4.356E-9,4.36E-9,4.364E-9,4.368E-9,4.372E-9,4.376E-9,4.38E-9,4.384E-9,4.388E-9,4.392E-9,4.396E-9,4.4E-9,4.404E-9,4.408E-9,4.412E-9,4.416E-9,4.42E-9,4.424E-9,4.428E-9,4.432E-9,4.436E-9,4.44E-9,4.444E-9,4.448E-9,4.452E-9,4.456E-9,4.46E-9,4.464E-9,4.468E-9,4.472E-9,4.476E-9,4.48E-9,4.484E-9,4.488E-9,4.492E-9,4.496E-9,4.5E-9,4.504E-9,4.508E-9,4.512E-9,4.516E-9,4.52E-9,4.524E-9,4.528E-9,4.532E-9,4.536E-9,4.54E-9,4.544E-9,4.548E-9,4.552E-9,4.556E-9,4.56E-9,4.564E-9,4.568E-9,4.572E-9,4.576E-9,4.58E-9,4.584E-9,4.588E-9,4.592E-9,4.596E-9,4.6E-9,4.604E-9,4.608E-9,4.612E-9,4.616E-9,4.62E-9,4.624E-9,4.628E-9,4.632E-9,4.636E-9,4.64E-9,4.644E-9,4.648E-9,4.652E-9,4.656E-9,4.66E-9,4.664E-9,4.668E-9,4.672E-9,4.676E-9,4.68E-9,4.684E-9,4.688E-9,4.692E-9,4.696E-9,4.7E-9,4.704E-9,4.708E-9,4.712E-9,4.716E-9,4.72E-9,4.724E-9,4.728E-9,4.732E-9,4.736E-9,4.74E-9,4.744E-9,4.748E-9,4.752E-9,4.756E-9,4.76E-9,4.764E-9,4.768E-9,4.772E-9,4.776E-9,4.78E-9,4.784E-9,4.788E-9,4.792E-9,4.796E-9,4.8E-9,4.804E-9,4.808E-9,4.812E-9,4.816E-9,4.82E-9,4.824E-9,4.828E-9,4.832E-9,4.836E-9,4.84E-9,4.844E-9,4.848E-9,4.852E-9,4.856E-9,4.86E-9,4.864E-9,4.868E-9,4.872E-9,4.876E-9,4.88E-9,4.884E-9,4.888E-9,4.892E-9,4.896E-9,4.9E-9,4.904E-9,4.908E-9,4.912E-9,4.916E-9,4.92E-9,4.924E-9,4.928E-9,4.932E-9,4.936E-9,4.94E-9,4.944E-9,4.948E-9,4.952E-9,4.956E-9,4.96E-9,4.964E-9,4.968E-9,4.972E-9,4.976E-9,4.98E-9,4.984E-9,4.988E-9,4.992E-9,4.996E-9,5E-9,5.004E-9,5.008E-9,5.012E-9,5.016E-9,5.02E-9,5.024E-9,5.028E-9,5.032E-9,5.036E-9,5.04E-9,5.044E-9,5.048E-9,5.052E-9,5.056E-9,5.06E-9,5.064E-9,5.068E-9,5.072E-9,5.076E-9,5.08E-9,5.084E-9,5.088E-9,5.092E-9,5.096E-9,5.1E-9,5.104E-9,5.108E-9,5.112E-9,5.116E-9,5.12E-9,5.124E-9,5.128E-9,5.132E-9,5.136E-9,5.14E-9,5.144E-9,5.148E-9,5.152E-9,5.156E-9,5.16E-9,5.164E-9,5.168E-9,5.172E-9,5.176E-9,5.18E-9,5.184E-9,5.188E-9,5.192E-9,5.196E-9,5.2E-9,5.204E-9,5.208E-9,5.212E-9,5.216E-9,5.22E-9,5.224E-9,5.228E-9,5.232E-9,5.236E-9,5.24E-9,5.244E-9,5.248E-9,5.252E-9,5.256E-9,5.26E-9,5.264E-9,5.268E-9,5.272E-9,5.276E-9,5.28E-9,5.284E-9,5.288E-9,5.292E-9,5.296E-9,5.3E-9,5.304E-9,5.308E-9,5.312E-9,5.316E-9,5.32E-9,5.324E-9,5.328E-9,5.332E-9,5.336E-9,5.34E-9,5.344E-9,5.348E-9,5.352E-9,5.356E-9,5.36E-9,5.364E-9,5.368E-9,5.372E-9,5.376E-9,5.38E-9,5.384E-9,5.388E-9,5.392E-9,5.396E-9,5.4E-9,5.404E-9,5.408E-9,5.412E-9,5.416E-9,5.42E-9,5.424E-9,5.428E-9,5.432E-9,5.436E-9,5.44E-9,5.444E-9,5.448E-9,5.452E-9,5.456E-9,5.46E-9,5.464E-9,5.468E-9,5.472E-9,5.476E-9,5.48E-9,5.484E-9,5.488E-9,5.492E-9,5.496E-9,5.5E-9,5.504E-9,5.508E-9,5.512E-9,5.516E-9,5.52E-9,5.524E-9,5.528E-9,5.532E-9,5.536E-9,5.54E-9,5.544E-9,5.548E-9,5.552E-9,5.556E-9,5.56E-9,5.564E-9,5.568E-9,5.572E-9,5.576E-9,5.58E-9,5.584E-9,5.588E-9,5.592E-9,5.596E-9,5.6E-9,5.604E-9,5.608E-9,5.612E-9,5.616E-9,5.62E-9,5.624E-9,5.628E-9,5.632E-9,5.636E-9,5.64E-9,5.644E-9,5.648E-9,5.652E-9,5.656E-9,5.66E-9,5.664E-9,5.668E-9,5.672E-9,5.676E-9,5.68E-9,5.684E-9,5.688E-9,5.692E-9,5.696E-9,5.7E-9,5.704E-9,5.708E-9,5.712E-9,5.716E-9,5.72E-9,5.724E-9,5.728E-9,5.732E-9,5.736E-9,5.74E-9,5.744E-9,5.748E-9,5.752E-9,5.756E-9,5.76E-9,5.764E-9,5.768E-9,5.772E-9,5.776E-9,5.78E-9,5.784E-9,5.788E-9,5.792E-9,5.796E-9,5.8E-9,5.804E-9,5.808E-9,5.812E-9,5.816E-9,5.82E-9,5.824E-9,5.828E-9,5.832E-9,5.836E-9,5.84E-9,5.844E-9,5.848E-9,5.852E-9,5.856E-9,5.86E-9,5.864E-9,5.868E-9,5.872E-9,5.876E-9,5.88E-9,5.884E-9,5.888E-9,5.892E-9,5.896E-9,5.9E-9,5.904E-9,5.908E-9,5.912E-9,5.916E-9,5.92E-9,5.924E-9,5.928E-9,5.932E-9,5.936E-9,5.94E-9,5.944E-9,5.948E-9,5.952E-9,5.956E-9,5.96E-9,5.964E-9,5.968E-9,5.972E-9,5.976E-9,5.98E-9,5.984E-9,5.988E-9,5.992E-9,5.996E-9,6E-9]

# To create a gaussian we need f = ae^(-((x-b)^2)/2c^2) - we need a, b, c.
# a = maximum height
# b = position of maximum
# c = FWHM / 2.35482

def findMaxAndHalfMax(x_points,y_points):
	base_level = float("inf")
	for y in y_points:
		if y < base_level:
			base_level = y
	#first we found the base_level at which to normalise all FWHMs
	a = 0
	a_index = 0
	for i, y in enumerate(y_points):
		if y > a:
			a_index = i
			a = y
			b = x_points[a_index]
	# We now have A and B; time to find c_half.
	a = a - base_level #to account for the base_level
	half_a = ((a)/2.0) + base_level #again, accounting for base
	delta = float("inf")
	c_half_index = 0
	# By finding the point with the smallest delta we find the point nearest the FWHM
	for i, y in enumerate(y_points[a_index:]):
		if abs(y-half_a) < delta:
			c_half_index = i+a_index
			delta = abs(y-half_a)
	c = ((x_points[c_half_index]-b)*2.0)/2.35482 #to convert from fwhm, hence 2.35482
	return base_level, a, b, c

def gaussian(base_level,a,b,c,x):
	return base_level + (a * np.exp(-(((x-b)**2)/(2*(c**2)))))

def generateGaussian(base_level, a, b, c, x_points):
	gaussian_set = []
	for x in x_points:
		gaussian_set.append(gaussian(base_level,a,b,c,x))
	return gaussian_set

def writeOut(data,filename):
	fileObject = open(filename,'w')
	csvObject = csv.writer(fileObject,delimiter=",")
	for row in data:
		csvObject.writerow(row)
	fileObject.close()

base_level, a, b, c = findMaxAndHalfMax(x_points,y_points)
gaussian_set = generateGaussian(base_level,a,b,c,x_points)

writeOut([gaussian_set],filename)
