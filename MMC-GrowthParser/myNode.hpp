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



#ifndef myNode_hpp
#define myNode_hpp
#include <vector>
#include <stdio.h>
#include <algorithm>


class myNode {
	double mTime;
	int mNumberOfDescendents;
	std::vector<int> mChildren;
    std::vector<int> mDescendents;

	int mLable;
public:
	//! empty constructor.
	myNode();
	//! destructor constructor.
	~myNode();
	//! A constructor seting mTime to T,  mLable to lable and the mNumberOfDescendents to NoD.
	/*!
	 \param t a double argument.
	 \param NoD an int argument
	 */
	myNode(double time, int lable, int NoD);
	//! A constructor seting mTime to T and mNumberOfDescendents and mLable to 0.
	/*!
	 \param t a double argument.
	 */
	myNode(double t, int lable);
	
	
	///////////////////////////////
	// adder, getter and setter //
	/////////////////////////////
    void addDescendent(int IndexDescendent);
    void sortDescendent();
	void addChild (myNode Node);
	std::vector<int> getChildren();
    std::vector<int> getDescendents();
	double getTime();
	int getNumberOfDescendents();
	int getLable();
	void setLable(int lable);
	void setTime(double time);
	
	
};



#endif /* myNode_hpp */
