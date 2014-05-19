#include "AdjacencyGraph.h"
#include "Genome.h"
#include "Chromosome.h"

using namespace std;

int main(int argc, char *argv[])
{
    Genome *a;
    Genome *b;
    
    try {
        a = new Genome("Genome A", "0691 0690 0689 0688 0687 -0686 -0685 -0684 0679 0678 -0675 -0674 -0673 0672 -0667 -0666 -0665 -0664 -0663 -0662 0661 -0660 -0659 -0658 -0655 1478 0653 -0652 0651 -0355 -0353 0351 -0348 -0347 0344 0343 0342 0341 0340 0337 -0335 -0334 0333 0332 0331 0330 -0329 -1477 0324 0323 0322 0321 0320 0319 0318 0317 0316 0315 -0313 0140 0139 0138 0137 -0119 0120 0121 0122 0123 -0124 -0125 -0126 0127 0128 0136 -1579 1118 1117 1116 -1114 -1113 -1112 1111 1110 1109 1107 1106 1105 1104 1103 1102 1101 1099 1098 1097 1096 1095 1094 1093 1090 1089 1087 1086 -1213 1212 1211 1210 1209 1207 1206 -1203 -1202 1201 1200 1199 1198 1197 1195 1194 1192 1191 1190 1189 1188 1187 1186 1185 1184 1183 1182 1181 1179 1178 1177 1176 -1174 -1173 -1172 -1171 -1170 -1169 -1168 -1167 -1166 -1165 -1164 -1163 1162 -1160 -1159 1158 1157 -1156 -1154 -1153 1152 -1151 -1150 1149 1148 1147 -1146 1144 -1142 -1141 -1140 -1139 -1138 -1137 -1136 -1135 -1132 -1131 1130 1127 1126 1125 1124 1123 -1122 1121 -0899 -0897 -0896 -0895 -0894 -0893 -0892 0891 -0890 0889 0887 0884 0882 0881 0880 -0878 -0875 -0874 0873 0872 0869 0868 0867 0866 0230 -0226 -0225 0224 0223 0222 0221 0220 0219 0218 0217 0216 0215 0214 -0406 -0405 0404 0403 0402 0401 0400 0399 -0396 -0395 0394 0393 0392 0391 0390 0389 -0387 -0386 -0384 -0383 0053 0382 0381 0380 0379 0378 0377 -0376 0374 0373 0372 0371 0370 0369 0368 0367 0365 0363 -0360 0357 -0156 -0059 -1001 -1002 -1003 -1006 -1007 -1008 -1009 -1010 1011 1012 1013 1014 1015 1016 1017 1022 1023 1024 -1025 -1027 1028 1032 1033 -1034 1035 -1036 1037 1038 1039 1043 1044 1045 1047 -1048 -1049 -1050 -1051 -1052 -1053 1054 1055 1058 1060 1061 0056 -1067 -1068 -1070 -1071 -1072 1073 1075 1076 1077 0060 1475 1474 1473 1080 1472 -1081 -1082 -1083 1084 -0113 0114 -0117 -0118 1234 1235 -1237 -1239 -1240 -1241 -1242 -1243 -1244 -1245 -1246 1248 1249 1250 -1251 1252 1254 1255 1434 -1256 1257 1258 1259 1260 1261 -1263 -1265 -1266 -1267 1271 1272 1273 1274 1275 1277 1278 1280 -1281 -1282 -1283 -1284 -1285 -1286 1287 1288 1289 1290 1291 1292 1293 -1294 -1296 -1297 -1299 1301 1302 1303 1304 1305 1306 1307 -1312 1313 1314 1315 1316 1317 1318 1319 -1321 -1322 -1323 -1324 -1325 -1326 1327 1328 1329 1330 -1332 -1333 -1336 -1337 -1338 1339 1340 0449 -1341 1342 1343 1345 1346 1347 1352 1353 1355 -1356 1357 0032 -0474 -0463 -0462 -0461 0458 0457 0455 0454 0453 -0452 -0451 0040 0188 0190 -0197 0198 0199 0201 -0202 -0203 -0204 -0207 1420 0141 0142 0143 -0144 -0145 -0147 -0148 -0149 -0150 -0151 -0153 -0518 0519 0521 -0525 0526 -0527 1386 -1385 -1384 1383 1382 1381 1380 1379 1378 1377 1376 1375 -1373 -1372 -1371 1364 1363 -1362 -1361 -0487 -0488 -0489 -0492 -0493 -0494 -0495 -0496 -0497 -0498 -0499 -0500 -0501 -0513 -0514 -0515 -0516 0590 -0589 -0588 -0586 0585 0584 0581 0579 0578 0577 0576 0575 -0573 -0572 -0571 -0568 -0567 0559 -0556 -0552 -0551 0550 0549 -0548 -0545 -0544 -0543 -0542 -0541 -0540 0539 -0534 -0533 -0296 -0294 0293 -0292 -0291 -0290 -0288 0287 0284 0037 0278 0277 1406 0276 0275 0274 0273 0272 -0830 -0836 -0837 -0838 -0839 -0840 0841 0842 0845 -0847 -0848 -0850 -0851 -0852 0853 0854 0855 0856 -0857 -0986 -0982 -0981 0980 0979 -0978 -0977 1404 0976 -0975 0169 0170 1575 0088 -0046 -0173 -0174 -0179 0180 0181 0182 -0735 -0736 -0737 -0738 -0739 -0740 -0741 -0742 -0743 -0744 -0745 -0746 -0747 -0748 -0035 -0749 -0750 -0751 -0752 -0753 -0754 -0755 -0756 -0757 -0758 -0759 -0760 -0761 -0762 0763 0765 0768 -0771 -0773 -0774 -0775 0776 0062 -0779 -0780 -0781 0782 0785 0787 -0791 -0792 -0796 -0798 -0799 -0800 0801 -0803 0804 0805 0806 0808 -0811 -0814 -0815 -0816 0819 0820 0821 0822 0823 0824 0825 0826 -0828 -0045 -0270 -0269 -0268 0412 -0012 -0436 -0447 -0286 -0285 -0004 -0010 1471 -1470 -0264 -0263 -0259 -0258 -0256 -0255 -0254 -0253 -0252 -0251 0250 0249 0248 0247 0246 0243 -0236 -0235 0234 -1216 -1217 -1218 -1219 -1220 -1221 1222 -1224 1225 -1226 1227 1228 1229 1230 1231 -0647 -0646 -0645 -0644 -0643 -0642 0641 0640 0635 0634 0633 0632 -0629 -0628 -0627 -0626 -0625 -0624 -0623 -0622 0621 0620 0619 0617 -0615 -0614 0613 0612 0611 0610 0609 0608 0607 0606 0605 0604 0603 0602 0601 0600 -0599 -0598 -0597 -0596 -0595 -0594 -0593 -0592 -0591 -0901 -0902 -0903 -0908 -0909 -0912 -0913 -0914 0915 0921 0922 -0923 -0924 0925 -0926 -0927 -0928 -0929 -0931 -1393 0935 0937 0938 -0939 -0940 -0941 -0945 -0946 -0947 0948 0950 -0951 0952 0953 0955 0956 0957 0959 -0960 -0961 -0963 0964 0966 -0967 -0968 -0969 -0970 -0971 -0972 -0973 -0974 0733 -0731 -0730 -0729 -0728 0727 -0726 -0725 0724 0723 0722 -0721 -0717 -0716 0715 0713 -0712 -0711 -0706 -0705 -0703 -0702 -0701 -0700 0697 0696 -0694 -0693 -0692 )");
        b = new Genome("Genome B", "0691 0690 0689 0688 0687 -0686 -0685 -0684 0679 0678 -0675 -0674 -0673 0672 -1206 -1207 -1209 -1210 -1211 -1212 1213 -1086 -1087 -1089 -1090 -1093 -1094 -1095 -1096 -1097 -1098 -1099 -1101 -1102 -1103 -1104 -1105 -1106 -1107 -1109 -1110 -1111 1112 1113 1114 -1116 -1117 -1118 1579 -0136 -0128 -0127 0126 0125 0124 -0123 -0122 -0121 -0120 0119 -0137 -0138 -0139 -0140 0313 -0315 -0316 -0317 -0318 -0319 -0320 -0321 -0322 -0323 -0324 0329 -0331 -0332 -0333 0334 0335 -0337 -0340 -0341 -0342 -0343 -0344 0347 0348 -0351 0353 0355 -0651 0652 0655 0659 0660 -0661 0662 0663 0664 0665 0666 0667 -1203 -1202 1201 1200 1199 1198 1197 1195 1194 1192 1191 1190 1189 1188 1187 1186 1185 1184 1183 1182 1181 1179 1178 1177 1176 -1174 -1173 -1172 -1171 -1170 -1169 -1168 -1167 -1166 -1165 -1164 -1163 1162 -1160 -1159 1158 1157 -1156 -1154 -1153 1152 -1151 -1150 1149 1148 1147 -1146 1144 -1142 -1141 -1140 -1139 -1138 -1137 -1136 -1135 -1132 -1131 1130 1127 1126 1125 1124 1123 -1122 1121 -0899 -0897 -0896 -0895 -0894 -0893 -0892 0891 -0890 0889 0887 0884 0881 0880 -0878 -0875 -0874 0873 0872 0869 0868 0867 0866 0230 -0226 -0225 0224 0223 0222 0221 0220 0219 0218 0217 0216 0215 0214 -0406 -0405 0404 0403 0402 0401 0400 0399 -0396 -0395 0394 0393 0392 0391 0390 0389 -0387 -0386 -0384 -0383 0053 0382 0381 0380 0379 0378 0377 -0376 0374 0373 0372 0371 0370 0369 0368 0367 0365 0363 -0360 0357 -0156 -0059 -1001 -1002 -1003 -1006 -1007 -1008 -1009 -1010 1011 1012 1013 1014 1015 1016 1017 1022 1023 1024 -1025 -1027 1028 1032 1033 -1034 1035 -1036 1037 1038 1039 1043 1044 1045 1047 -1048 -1049 -1050 -1051 -1052 -1053 1054 1055 1058 1060 1061 0056 -1067 -1068 -1070 -1071 -1072 1073 1075 1076 1077 0060 1475 1474 1473 1080 1472 -1081 -1082 -1083 1084 -0113 0114 -0117 -0118 1234 1235 -1237 -1239 -1240 -1241 -1242 -1243 -1244 -1245 -1246 1248 1249 1250 -1251 1252 1254 1255 1434 -1256 1257 1258 1259 1260 1261 -1263 -1265 -1266 -1267 1271 1272 1273 1274 1275 1277 1278 1280 -1281 -1282 -1283 -1284 -1285 -1286 1287 1288 1289 1290 1291 1292 1293 -1294 -1296 -1297 -1299 1301 1302 1303 1304 1305 1306 1307 -1312 1313 1314 1315 1316 1317 1318 1319 -1321 -1322 -1323 -1324 -1325 -1326 1327 1328 1329 1330 -1332 -1333 -1336 -1337 -1338 1339 1340 0449 -1341 1342 1343 1345 1346 1347 1352 1353 1355 -1356 1357 0032 -0474 -0463 -0462 -0461 0458 0457 0455 0454 0453 -0452 -0451 0040 0188 0190 -0197 0198 0199 0201 -0202 -0203 -0204 -0207 1420 0141 0142 0143 -0144 -0145 -0147 -0148 -0149 -0150 -0153 -0518 0519 0521 -0525 0526 -0527 1386 -1385 -1384 1383 1382 1381 1380 1379 1378 1377 1376 1375 -1373 -1372 -1371 1364 1363 -1362 -1361 -0487 -0488 -0489 -0492 -0493 -0494 -0495 -0496 -0497 -0498 -0499 -0500 -0501 -0513 -0514 -0515 -0516 0590 -0589 -0588 -0586 0585 0584 0581 0579 0578 0577 0576 0575 -0573 -0572 -0571 -0568 -0567 0559 -0556 -0552 -0551 0550 0549 -0548 -0545 -0544 -0543 -0542 -0541 -0540 0539 -0534 -0533 -0296 -0294 0293 -0292 -0291 -0290 -0288 0287 0284 0037 0278 0277 1406 0276 0275 0274 0273 0272 -0830 -0836 -0837 -0838 -0839 -0840 0841 0842 0845 -0847 -0848 -0850 -0851 -0852 0853 0854 0855 0856 -0857 -0986 -0982 -0981 0980 0979 -0978 -0977 1404 0976 -0975 0169 0170 0088 -0046 -0173 -0174 -0179 0180 0181 0182 -0184 -0735 -0736 -0737 -0738 -0739 -0740 -0741 -0742 -0743 -0744 -0745 -0746 -0747 -0748 -0035 -0749 -0750 -0751 -0752 -0753 -0754 -0755 -0756 -0757 -0758 -0759 -0760 -0761 -0762 0763 0765 0768 -0771 -0773 -0774 -0775 0776 0062 -0779 -0780 -0781 0782 0785 0787 -0791 -0792 -0796 -0798 -0799 -0800 0801 -0803 0804 0805 0808 -0811 -0814 -0815 -0816 0819 0820 0821 0822 0823 0824 0825 0826 -0828 -0045 -0270 -0269 -0268 -0264 -0263 -0259 -0258 -0256 -0255 -0254 -0253 -0252 -0251 0250 0249 0248 1519 0247 0246 0243 -0236 -0235 0234 -1216 -1217 -1218 -1219 -1220 -1221 1222 -1224 1225 -1226 1227 1228 1229 1230 1231 -0647 -0646 -0645 -0644 -0643 -0642 0641 0640 0635 0634 0633 0632 -0629 -0628 -0627 -0626 -0625 -0624 -0623 -0622 0621 0620 0619 0617 -0615 -0614 0613 0612 0611 0610 0609 0608 0607 0606 0605 0604 0603 0602 0601 0600 -0599 -0598 -0597 -0596 -0595 -0594 -0593 -0592 -0591 -0901 -0902 -0903 -0908 -0909 -0912 -0913 -0914 0915 0921 0922 -0923 -0924 0925 -0926 -0927 -0928 -0929 -0931 0935 0937 0938 -0939 -0940 -0941 -0945 -0946 -0947 0948 0950 -0951 0952 0953 0955 0956 0957 0959 -0960 -0961 -0963 0964 0966 -0967 -0968 -0969 -0970 -0971 -0972 -0973 -0974 0733 -0731 -0730 -0729 -0728 0727 -0726 -0725 0724 0723 0722 -0721 -0717 -0716 0715 0713 -0712 -0711 -0706 -0705 -0703 -0702 -0701 -0700 0697 0696 -0694 -0693 -0692 )");
        //a = new Genome("Genome A", "1 19 5 -20 3 | 4 21 2 ) 26 -27 )");
        //b = new Genome("Genome B", "1 23 24 2 3 25 4 5 |");
        //a = new Genome("Genome A", "1 19 5 -20 3 | 4 21 2 ) 26 -27 )");
        //b = new Genome("Genome B", "1 23 24 2 3 25 4 5 |");
    } catch (const char *err) {
        std::cerr << err << std::endl;
        return 1;
    }
/*
    		    Chromosome *a1 = new Chromosome("chrA1", true);

		    a1->genes.push_back(-1);
		    a1->genes.push_back(19);

		    Chromosome *a2 = new Chromosome("chrA2", true);

		    a2->genes.push_back(2);
		    a2->genes.push_back(20);
		    a2->genes.push_back(3);

                    Chromosome *a3 = new Chromosome("chrA3", true);

                    a3->genes.push_back(-7);
                    a3->genes.push_back(-4);
                    a3->genes.push_back(6);

                    Chromosome *a4 = new Chromosome("chrA4", true);

                    a4->genes.push_back(5);
                    a4->genes.push_back(23);
                    a4->genes.push_back(24);
                    a4->genes.push_back(9);
                    a4->genes.push_back(-11);
                    a4->genes.push_back(12);

                    Chromosome *a5 = new Chromosome("chrA5", true);

                    a5->genes.push_back(-8);
                    a5->genes.push_back(-10);
                    a5->genes.push_back(28);

		    Genome *a = new Genome("Genome A");
		    a->chromosomes.push_back(a1);
		    a->chromosomes.push_back(a2);
                    a->chromosomes.push_back(a3);
                    a->chromosomes.push_back(a4);
                    a->chromosomes.push_back(a5);

                    Chromosome *b1 = new Chromosome("chrB1", true);

		    b1->genes.push_back(-1);
		    b1->genes.push_back(-2);
		    b1->genes.push_back(-3);
		    b1->genes.push_back(21);
		    b1->genes.push_back(-4);
		    b1->genes.push_back(-7);
		    b1->genes.push_back(22);
		    b1->genes.push_back(-6);

                    Chromosome *b2 = new Chromosome("chrB2", true);

		    b2->genes.push_back(-5);
		    b2->genes.push_back(-9);
		    b2->genes.push_back(25);
		    b2->genes.push_back(11);
		    b2->genes.push_back(-12);
		    b2->genes.push_back(26);
		    b2->genes.push_back(8);
		    b2->genes.push_back(10);
                    b2->genes.push_back(27);

		    Genome *b = new Genome("Genome B");
		    b->chromosomes.push_back(b1);
                    b->chromosomes.push_back(b2);
*/
		    AdjacencyGraph *ag = new AdjacencyGraph(a,b);

                    std::queue<Genome *> steps;
                    std::queue< Rearrangement > dcjs;
		    
                    std::cout << "Distancia pelo SortingByDCJ: " << ag->sortByDCJsubst(steps, dcjs) << std::endl;

                    //std::cout << "Distancia pela fórmula: " << ag->DCJdistance() << std::endl;

                    //std::cout << "Distancia pelo SortingByDCJRestrict: "
                            //<< ag->sortByRestrictedDCJ() << std::endl;
}
