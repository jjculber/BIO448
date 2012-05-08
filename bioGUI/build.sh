#!/bin/bash
make clean
make
jar cvfm BioProgram.jar manifest.txt bioGUI/*
