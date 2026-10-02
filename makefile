CXX = g++
CXXFLAGS = -g -std=c++17 -Wall -Wextra -Iinclude $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --glibs)

OBJDIR = objs
BINDIR = bin
SRCDIR = src

.DEFAULT_GOAL := all
.PHONY: all calibration clean

OBJS = \
  $(OBJDIR)/main.o \
  $(OBJDIR)/LineShapes.o \
  $(OBJDIR)/FitConfig.o \
  $(OBJDIR)/FitModel.o \
  $(OBJDIR)/FitResults.o \
  $(OBJDIR)/SpectrumIO.o \
  $(OBJDIR)/Plotting.o \
  $(OBJDIR)/UpperLimit.o \
  $(OBJDIR)/SpectrumCalibration.o

CALIBRATION_SOURCES = $(SRCDIR)/SpectrumCalibration.cpp $(SRCDIR)/LineShapes.cpp
CALIBRATION_HEADERS = include/SpectrumCalibration.h include/LineShapes.h

# A plain `make` builds both the existing fitter and the calibration command.
all: $(BINDIR)/bwgaussfit $(BINDIR)/bwcalibrate

# Build just the standalone calibration command with `make calibration`.
calibration: $(BINDIR)/bwcalibrate

$(BINDIR)/bwgaussfit: $(OBJS)
	@mkdir -p $(BINDIR)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(ROOTLIBS)

# This define enables the standalone main() only for bwcalibrate.
$(BINDIR)/bwcalibrate: $(CALIBRATION_SOURCES) $(CALIBRATION_HEADERS)
	@mkdir -p $(BINDIR)
	$(CXX) $(CXXFLAGS) -DBWFIT_CALIBRATION_MAIN $(CALIBRATION_SOURCES) -o $@ $(ROOTLIBS)

$(OBJDIR)/SpectrumCalibration.o: $(CALIBRATION_HEADERS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJDIR)/*.o $(BINDIR)/bwgaussfit $(BINDIR)/bwcalibrate
