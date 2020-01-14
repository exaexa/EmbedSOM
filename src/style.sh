#!/bin/sh
clang-format -style="{BasedOnStyle: Mozilla, UseTab: ForIndentation, IndentWidth: 8, TabWidth: 8, AccessModifierOffset: -8, PointerAlignment: Left}" -verbose -i *.cpp *.h
