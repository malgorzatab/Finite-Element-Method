#include "Node.h"

Node::Node() {
	this->x = 0;
	this->y = 0;
	this->temperature = 0;
}

Node::Node(double x, double y, double temperature, bool status) {
	this->x = x;
	this->y = y;
	this->temperature = temperature;
	this->status = status;

}

Node::~Node(){}


int Node::getId() {
	return this->id;
}

double Node::getX() {
	return this->x;
}

double Node::getY() {
	return this->y;
}

double Node::getTemperature() {
	return this->temperature;
}

bool Node::getStatus() {
	return this->status;
}


void Node::setTemperature(double temperature) {
	this->temperature = temperature;
}

void Node::setStatus(bool status) {
	this->status = status;
}

void Node::setX(double x) {
	this->x = x;
}

void Node::setY(double y) {
	this->y = y;
}

void Node::setId(int id) {
	this->id = id;
}