cd ./lm/
echo "*************************linear switcher***************************************"
echo "*************************linear switcher STC***************************************"
spaceex -g 5_dim_linear_switch.cfg -m 5_dim_linear_switch.xml -v l -o out.gen

echo "*************************linear switcher LGG***************************************"
spaceex -g lmlgg.cfg -m 5_dim_linear_switch.xml -v l -o out.gen
cd ../
cd ./platoon/
echo "*************************Platoon***************************************"
echo "*************************PLAD01-BND42STC***************************************"
spaceex -g PLAD01-BNDSTC.cfg -m PLAD01-BND.xml -v l -o out.gen
echo "*************************PLAD01-BND42LGG***************************************"
spaceex -g PLAD01-BND42.cfg -m PLAD01-BND.xml -v l -o out.gen
echo "*************************PLAN01-UNBSTC***************************************"
spaceex -g PLAN01-UNBSTC.cfg -m PLAN01-UNB.xml -v l -o out.gen
echo "*************************PLAD01-UNBLGG***************************************"
spaceex -g PLAN01-UNB50.cfg -m PLAN01-UNB.xml -v l -o out.gen
cd ../
cd ./spacecraft/
echo "*************************SpaceCraft***************************************"
echo "************************* SRA01-SR0 STC***************************************"
spaceex -g SRA01-SRSTC.cfg -m SRA01-SR0_.xml -v l -o out.gen
echo "************************* SRA01-SR0 LGG***************************************"
spaceex -g SRA01-plot.cfg -m SRA01-SR0_.xml -v l -o out.gen
echo "************************* SRNA01-SR0 STC***************************************"
spaceex -g SRNA01-SR0STC.cfg -m SRNA01-SR0_.xml -v l -o out.gen
echo "************************* SRNA01-SR0 LGG***************************************"
spaceex -g SRNA01-SR01.cfg -m SRNA01-SR0_.xml -v l -o out.gen
cd ../
