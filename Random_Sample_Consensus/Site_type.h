#ifndef SITE_TYPE_H
#define SITE_TYPE_H

class Site
{
public:
	void set_x(double x) { this->x = x; };
	void set_y(double y) { this->y = y; };
	double get_x() { return this->x; };
	double get_y() { return this->y; };
private:
	double x;
	double y;
}; 

#endif // !SITE_TYPE_H

