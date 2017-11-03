/*
 * MMC_Growth is used to simulate gene trees under the psi- and kingman
 * coalescent process with exponential growth.
 *
 * Copyright (C) 2016 Sebastian Matuszewski & Marcel Hildebrandt
 *
 * This file is part of MMC_Growth.
 *
 * MMC_Growth is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef mySummaryStatistics_hpp
#define mySummaryStatistics_hpp

#include <stdio.h>
#include "myTree.hpp"

class mySummaryStatistics{
    double mThetaW;     // Watterson's Theta
    double mThetaT;     // Tajima's Theta
    double mThetaH;     // Fay and Wu's Theta
    double mTajimasD;   // Tajima's D
    double mFuAndLisD;  // Fu and Li's D
    double mFayAndWusH; // Fay and Wu's H
    
    public:
    
    //  constructor
    mySummaryStatistics(myTree tree);
    
    
    /////////////
    // getter //
    ///////////
    double getTheta_W();
    double getThetaT();
    double getTajimasD();
    double getFuAndLisD();
    double getThetaH();
    double getFayAndWusH();
    
};

template <typename _typeName>
std::string convert_NumberToStringTree ( _typeName Number );


#endif /* mySummaryStatistics_hpp */
