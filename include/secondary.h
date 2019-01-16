#ifndef INCLUDE_SECONDARY_H_
#define INCLUDE_SECONDARY_H_

class SecondarySource {
public:
	SecondarySource();
	virtual ~SecondarySource();
	double get(const double& T) const;

protected:
	int A = 0;
	int Z = 0;
};


#endif /* INCLUDE_SECONDARY_H_ */
