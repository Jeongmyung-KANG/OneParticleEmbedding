
#g++ MultiMixing.cxx -o m -I/Users/gangjeongmyeong/RooUnfold-1.1.1/src `fastjet-config --cxxflags --libs --plugins` `root-config --cflags` `root-config --glibs` -L/Users/gangjeongmyeong/RooUnfold-1.1.1 -lRooUnfold  

#g++ Mixing.cxx -o sj -I/Users/gangjeongmyeong/RooUnfold-1.1.1/src `fastjet-config --cxxflags --libs --plugins` `root-config --cflags` `root-config --glibs` -L/Users/gangjeongmyeong/RooUnfold-1.1.1 -lRooUnfold  

#g++ Save_SE.cc -o sr -I/Users/gangjeongmyeong/RooUnfold-1.1.1/src `fastjet-config --cxxflags --libs --plugins` `root-config --cflags` `root-config --glibs` -L/Users/gangjeongmyeong/RooUnfold-1.1.1 -lRooUnfold  
#g++ Save_ME.cc -o ME -I/Users/gangjeongmyeong/RooUnfold-1.1.1/src `fastjet-config --cxxflags --libs --plugins` `root-config --cflags` `root-config --glibs` -L/Users/gangjeongmyeong/RooUnfold-1.1.1 -lRooUnfold  



#g++ Save_SE.cc -o sr -I/Users/gangjeongmyeong/RooUnfold-1.1.1/src -I/Users/gangjeongmyeong/fjcontrib-1.100/ConstituentSubtractor `fastjet-config --cxxflags --libs --plugins` `root-config --cflags` `root-config --glibs` -L/Users/gangjeongmyeong/RooUnfold-1.1.1 -lRooUnfold -L/Users/gangjeongmyeong/fjcontrib-1.100/ConstituentSubtractor -LConstituentSubtractor
#g++ Save_ME.cc -o ME -I/Users/gangjeongmyeong/RooUnfold-1.1.1/src -I/Users/gangjeongmyeong/fjcontrib-1.100/ConstituentSubtractor `fastjet-config --cxxflags --libs --plugins` `root-config --cflags` `root-config --glibs` -L/Users/gangjeongmyeong/RooUnfold-1.1.1 -lRooUnfold -L/Users/gangjeongmyeong/fjcontrib-1.100/ConstituentSubtractor -LConstituentSubtractor
g++ SaveEvents.cxx -o sr -I/opt/homebrew/Cellar/fjcontrib/1.055/include -I/Users/gangjeongmyeong/RooUnfold-1.1.1/src `fastjet-config --cxxflags --libs --plugins` `root-config --cflags` `root-config --glibs` -L/Users/gangjeongmyeong/RooUnfold-1.1.1 -lRooUnfold -L/opt/homebrew/Cellar/fjcontrib/1.055/lib -lConstituentSubtractor
