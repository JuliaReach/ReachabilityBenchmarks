cd /Users/kpotomkin/Documents/MATLAB/SpaceEXHSCC2020/fo
cd ./256/
#echo "*************************256***************************************"
#spaceex -g filtered_oscillator.256.cfg -m filtered_oscillator_256.xml -v l -o out.gen
#echo "*************************256stc***************************************"
#spaceex -g filtered_oscillator.256stc.cfg -m filtered_oscillator_256.xml -v l -o out.gen

echo "*************************512***************************************"
spaceex -g filtered_oscillator.512.cfg -m filtered_oscillator_256.xml -v l -o out.gen
echo "*************************512stc***************************************"
spaceex -g filtered_oscillator.512stc.cfg -m filtered_oscillator_256.xml -v l -o out.gen
cd ../


