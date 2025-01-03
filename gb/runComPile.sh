
#g++ MultiMixing.cxx -o m -I/Users/gangjeongmyeong/RooUnfold-1.1.1/src `fastjet-config --cxxflags --libs --plugins` `root-config --cflags` `root-config --glibs` -L/Users/gangjeongmyeong/RooUnfold-1.1.1 -lRooUnfold  

#g++ Mixing.cxx -o sj -I/Users/gangjeongmyeong/RooUnfold-1.1.1/src `fastjet-config --cxxflags --libs --plugins` `root-config --cflags` `root-config --glibs` -L/Users/gangjeongmyeong/RooUnfold-1.1.1 -lRooUnfold  
//g++ Save_SE.cc -o sr -I/Users/gangjeongmyeong/RooUnfold-1.1.1/src `fastjet-config --cxxflags --libs --plugins` `root-config --cflags` `root-config --glibs` -L/Users/gangjeongmyeong/RooUnfold-1.1.1 -lRooUnfold  

 
g++ Save_ME.cc /src/GenEvents.cxx-o ME -I/Users/gangjeongmyeong/RooUnfold-1.1.1/src -I/include `fastjet-config --cxxflags --libs --plugins` `root-config --cflags` `root-config --glibs` -L/Users/gangjeongmyeong/RooUnfold-1.1.1 -lRooUnfold  
