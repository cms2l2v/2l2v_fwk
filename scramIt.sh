#!/bin/bash

USER_CXXFLAGS="-Wno-error=delete-non-virtual-dtor -Wno-error=unused-but-set-variable -Wno-error=unused-variable -Wno-error=sign-compare -Wno-error=reorder" scram b -j 20