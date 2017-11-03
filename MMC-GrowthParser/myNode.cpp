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



#include "myNode.hpp"

myNode::myNode(double time, int lable, int NoD)
{
	mTime = time;
	mLable = lable;
	mNumberOfDescendents = NoD;
    
}

myNode::myNode(double time, int lable)
{
	mLable = lable;
	mTime = time;
	mNumberOfDescendents = 0;
}


myNode::~myNode()
{}


myNode::myNode()
:mNumberOfDescendents()
{
	mTime = -1;
}


void myNode::addChild (myNode Node)
{
	this->mChildren.push_back(Node.getLable());
	this->mNumberOfDescendents += Node.mNumberOfDescendents;
    for (int i = 0; i < Node.mNumberOfDescendents; ++i) {
        this->mDescendents.push_back(Node.mDescendents[i]);
    }
}

void myNode::sortDescendent()
{
    std::sort(this->mDescendents.begin(), this->mDescendents.end());
}

std::vector<int> myNode::getChildren()
{
	return this->mChildren;
}


void myNode::addDescendent(int IndexDescendent)
{
    this->mDescendents.push_back(IndexDescendent);
}


std::vector<int> myNode::getDescendents()
{
    return this->mDescendents;
}


double myNode::getTime()
{
	return this->mTime;
}


int myNode::getNumberOfDescendents()
{
	return this->mNumberOfDescendents;
}


int myNode::getLable()
{
	return this->mLable;
}


void myNode::setLable(int lable)
{
	mLable = lable;
}


void myNode::setTime(double time)
{
	this->mTime = time;
}



