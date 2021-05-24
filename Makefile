CXX = g++
 
CXXFLAGS  = -Wall -std=c++14
 
TARGET = SpheCow

SRCS = BPLModel.cpp DeVaucouleursModel.cpp DensityModel.cpp EinastoModel.cpp GammaModel.cpp GaussLegendre.cpp HernquistModel.cpp HypervirialModel.cpp IsochroneModel.cpp JaffeModel.cpp Model.cpp NFWModel.cpp NukerModel.cpp PlummerModel.cpp PerfectSphereModel.cpp SersicModel.cpp SigmoidDensityModel.cpp SigmoidSurfaceDensityModel.cpp SpheCow.cpp SurfaceDensityModel.cpp ZhaoModel.cpp

OBJS=$(subst .cpp,.o,$(SRCS))
 
all: $(TARGET)
 
$(TARGET): $(OBJS)
	$(CXX) $(CXXLAGS) -o $(TARGET) $(OBJS)
 
clean:
	$(RM) $(TARGET)