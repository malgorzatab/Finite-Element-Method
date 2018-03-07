#pragma once
class Node {

public:

	int id;
	double x, y, temperature;		//ta temp = t0 z grid w ka¿dym wezle
	bool status; 

	Node();
	Node(double x, double y, double temperature, bool status);

	~Node();

	int getId();
	double getX();
	double getY();
	double getTemperature();
	bool getStatus();
	
	void setTemperature(double temperature);
	void setStatus(bool status);
	void setX(double x);
	void setY(double y);
	void setId(int id);


};