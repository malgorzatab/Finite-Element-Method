#pragma once
class Element {

public:
	int id_element;
	int id_nodes[4];
	int id_walls[4];

	double** lH;		//4x4

	double lP[4];

	void set_localH(double**lH);
	void set_localP(double*lP);

	Element();
	Element(int,int*,int*,double**, double*);
	~Element();

	int getId_element();
	int getId_nodes();
	int getId_walls();

	void setId_element(int id_element);
	void setId_nodes(int* id_nodes);
	void setId_walls(int* id_walls);


	void setlH(double**lH);
	void setlP(double*lP);


};
