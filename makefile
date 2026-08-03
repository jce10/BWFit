CXX = g++
CXXFLAGS = -g -std=c++17 -Wall -Wextra -Iinclude $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --glibs)

OBJDIR = objs
BINDIR = bin
SRCDIR = src

OBJS = \
  $(OBJDIR)/main.o \
  $(OBJDIR)/LineShapes.o \
  $(OBJDIR)/FitConfig.o \
  $(OBJDIR)/FitModel.o \
  $(OBJDIR)/FitResults.o \
  $(OBJDIR)/SpectrumIO.o \
  $(OBJDIR)/Plotting.o \
  $(OBJDIR)/UpperLimit.o

$(BINDIR)/bwgaussfit: $(OBJS)
	@mkdir -p $(BINDIR)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(ROOTLIBS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJDIR)/*.o $(BINDIR)/bwgaussfit
