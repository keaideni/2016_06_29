#include "OP.h"

struct Eigstruct
{
	int q;
	double lamda;
	VectorXd state;
};


std::string itos(int i)
{
        std::stringstream s;
        s << i;
        return s.str();
}

bool comp(const Eigstruct& a, const Eigstruct& b);
bool comp(const Eigstruct& a, const Eigstruct& b)
{
	return (a.lamda > b.lamda);
}

OP::OP(){};

OP::~OP(){};

OP::OP(const int& Dmin, const int& Dmax, const int& Qmin, const int& Qmax, const int& delta, const int& str)
{
	iniRaiseQ(Dmin, Dmax, Qmin, Qmax, delta, str);
}

OP::OP(const OP& a)
{
	QDim = a.QDim;
	QMat = a.QMat;
	RLQ = a.RLQ;
}


//delta for innihilation is -1, and for creation is 1. the N and I operators' delta is 0.
void OP::iniRaiseQ(const int& Dmin, const int& Dmax, const int& Qmin, const int& Qmax, const int& delta, const int& str)
{
	for (int i = Dmin; i <= Dmax; i++)
	{
		
		QDim.insert(std::pair<int, int>(i, 1));
	}

	for (int i = Qmin; i <= Qmax; i++)
	{
		
		RLQ.insert(std::pair<int, int>(i, i + delta));
		double y(0.0);

		getValue(str, i, y);

		MatrixXd tempm(1, 1);
		tempm << y;

		QMat.insert(std::pair<int, MatrixXd>(i, tempm));
	}
}


void OP::iniRaiseQ(const int str)
{


	

	MatrixXd temp(1, 1);
	

	QDim.insert(std::pair<int, int>(0, 1));
	QDim.insert(std::pair<int, int>(1, 1));


	switch (str)
	{

		case 1:
		{
		      RLQ.insert(std::pair<int, int>(0, 0));
		      RLQ.insert(std::pair<int, int>(1, 1));

		      temp(0, 0) = -1;

		      QMat.insert(std::pair<int, MatrixXd>(0, temp));


		      temp(0, 0) = 1;
		      QMat.insert(std::pair<int, MatrixXd>(1, temp));


		      break;
		}
		case 2:
		{
		      RLQ.insert(std::pair<int, int>(0, 1));

		      temp(0, 0) = 1;

		      QMat.insert(std::pair<int, MatrixXd>(0, temp));

		      break;
		}
		case 3:
		{
		      RLQ.insert(std::pair<int, int>(1, 0));

		      temp(0, 0) = 1;

		      QMat.insert(std::pair<int, MatrixXd>(1, temp));


		      break;
		}
		case 4:
		{
		      RLQ.insert(std::pair<int, int>(0, 0));
		      RLQ.insert(std::pair<int, int>(1, 1));

		      temp(0, 0) = 1;

		      QMat.insert(std::pair<int, MatrixXd>(0, temp));
		      QMat.insert(std::pair<int, MatrixXd>(1, temp));
		      break;
		}
	}
}







void OP::getValue(const int& str, const int& i, double& y)
{
	switch (str)
	{
		//for boson annihilation.
		case 1:
		{
		      y = sqrt((double)i);
		      break;
		}
		//for boson creation
		case 2:
		{
		      y = sqrt((double)(i + 1));
		      break;
		}
		//for I.
		case 3:
		{
		      y = 1.0;
		      break;
		}
		//for n.
		case 4:
		{
		      y = (double)i;
		      break;
		}
	}
}




void OP::findDim(const OP& a, const OP& b, std::unordered_map<int, int> &oldDim, std::unordered_map<std::pair<int, int>, int, classcom> &startDim) const
{
	for (int na = 0; na <= OP::Max; ++na)
	{
                auto ita = a.QDim.find(na);
                if(ita == a.QDim.end()) continue;
		for (int nb = 0; nb <= OP::Max; ++nb)
		{
                        auto itb = b.QDim.find(nb);
                        if(itb == b.QDim.end()) continue;
			int tempQ;
			tempQ=ita->first + itb->first;

			bool i = (tempQ <= Max);
			//bool j = ((a.RLQ.at(ita->first).QN + b.RLQ.at(itb->first).QN) <= Max);

			if (i)
			{
				auto  itc = oldDim.find(tempQ);
				if (itc == oldDim.end())
				{
					oldDim.insert(std::pair<int, int>(tempQ, ita->second * itb->second));
					startDim[std::pair<int, int>(ita->first, itb->first)] = 0;
				}
				else
				{
					startDim[std::pair<int, int>(ita->first, itb->first)] = oldDim[itc->first];
					oldDim[itc->first] += ita->second * itb->second;
				}
			}
		}

	}
}

void OP::kronO(const OP &a, const OP &b)
{
	clear();





	//find the dimention of each good quantum number, and label the place to put the kron of two blocks.
	std::unordered_map<std::pair<int, int>, int, classcom> startDim;
	findDim(a, b, QDim, startDim);



	//first calculate the good quantum number.
	for (auto ita = a.RLQ.begin(); ita != a.RLQ.end(); ita++)
	{
		for (auto itb = b.RLQ.begin(); itb != b.RLQ.end(); itb++)
		{
			int tempQR;
			tempQR = ita->first + itb->first;

			bool i = (tempQR <= Max);
			bool j = ((a.RLQ.at(ita->first) + b.RLQ.at(itb->first)) <= Max);

			if (i&&j)
			{



				auto itc = RLQ.find(tempQR);
				if (itc == RLQ.end())
				{
					int tempQL;
					tempQL = ita->second + itb->second;
					RLQ.insert(std::pair<int, int>(tempQR, tempQL));



					int col = QDim.at(tempQR);
					int row = QDim.at(tempQL);

					//creat a zero mat matched the dimention of good quantum number.
					MatrixXd tempm(MatrixXd::Zero(row, col));
					
					QMat.insert(std::pair<int, MatrixXd>(tempQR, tempm));
				}

			}
		}
	}









	//calculate the kron and put it in the place.
	for (auto ita = a.QMat.begin(); ita != a.QMat.end(); ita++)
	{
		for (auto itb = b.QMat.begin(); itb != b.QMat.end(); itb++)
		{




			int tempQR;
			tempQR = ita->first + itb->first;

			bool i = (tempQR <= Max);
			bool j = ((a.RLQ.at(ita->first) + b.RLQ.at(itb->first) )<= Max);

			if (i&&j)
			{


				

				


				int startR = startDim[std::pair<int, int>(ita->first, itb->first)];




				int startL = startDim[std::pair<int, int>(a.RLQ.at(ita->first), b.RLQ.at(itb->first))];



				for (int m = 0; m < ita->second.rows(); m++)
				{
					for (int n = 0; n < ita->second.cols(); n++)
					{
						int a = itb->second.rows();
						int b = itb->second.cols();
						
						//QMat.at(tempQR).block(startL + m*itb->second.rows(), startR + n*itb->second.cols(),a,b)
						//	= ita->second(m, n)*itb->second;
						for (int i = 0; i < itb->second.rows(); ++i)
						{
							for (int j = 0; j < itb->second.cols(); ++j)
							{
								QMat.at(tempQR)(startL + m*itb->second.rows() + i, startR + n*itb->second.cols() + j)
										= ita->second(m, n)*itb->second(i,j);
							}
						}
						
					}
				}



				
			}

		}
	}
}



void OP::transO(const OP& a)
{
	QDim = a.QDim;

	for (auto it = a.RLQ.begin(); it != a.RLQ.end(); it++)
	{
		RLQ.insert(std::pair<int,int>(it->second, it->first));
		QMat.insert(std::pair<int, MatrixXd>(it->second, a.QMat.at(it->first).transpose()));
	}

}


void OP::add(const OP&a, const OP& b)
{
	if (a.QDim.size() != b.QDim.size())
	{
		std::cout << "the Q dim size doesn't match!" << std::endl;
		exit(1);
	}

	for (auto tempit = a.QDim.begin(); tempit != a.QDim.end(); tempit++)
	{
		if (tempit->second != b.QDim.at(tempit->first))
		{
			std::cout << "the Q dimention doesn't match each other!" << std::endl;
			exit(1);
		}
	}

	QDim = a.QDim;
	RLQ = a.RLQ;
	QMat = a.QMat;





	for (auto tempit = b.QMat.begin(); tempit != b.QMat.end(); tempit++)
	{
		auto it = QMat.find(tempit->first);
		if (it == QMat.end())
		{
			QMat.insert(std::pair<int, MatrixXd>(tempit->first, b.QMat.at(tempit->first)));
			RLQ.insert(std::pair<int, int>(tempit->first, b.RLQ.at(tempit->first)));

		}
		else
		{
			if (a.RLQ.at(tempit->first) != b.RLQ.at(tempit->first))
			{
				std::cout << "add0: RQ != LQ, WRONG!" << std::endl;
				exit(1);
			}
			QMat.at(tempit->first) += b.QMat.at(tempit->first);
		}
	}



}






void OP::add(const OP& a)
{
	if (QDim.size() != a.QDim.size())
	{
		std::cout << "the Q dim size doesn't match!" << std::endl;
		exit(1);
	}

	for (auto tempit = QDim.begin(); tempit != QDim.end(); tempit++)
	{
		if (tempit->second != a.QDim.at(tempit->first))
		{
			std::cout << "the Q dimention doesn't match each other!" << std::endl;
			exit(1);
		}
	}
	for (auto tempit = a.QMat.begin(); tempit != a.QMat.end(); tempit++)
	{
		auto it = QMat.find(tempit->first);
		if (it == QMat.end())
		{
			QMat.insert(std::pair<int, MatrixXd>(tempit->first, a.QMat.at(tempit->first)));
			RLQ.insert(std::pair<int, int>(tempit->first, a.RLQ.at(tempit->first)));
		}

		else
		{
			/*if(RLQ.at(tempit->first).QN != a.RLQ.at(tempit->first).QN)
			{
			std::cout<<"add1: RQ != LQ, WRONG!"<<std::endl;
			exit(1);
			}*/

			QMat.at(tempit->first) += a.QMat.at(tempit->first);

		}
	}
}



//reride the operator + to calculate the sum of two operator.
OP OP::operator +(const OP& a)
{
	OP sum;


	//test the QDim.
	if (QDim.size() != a.QDim.size())
	{
		std::cout << "the Q dim size doesn't match!" << std::endl;
		exit(1);
	}
	for (auto tempit = QDim.begin(); tempit != QDim.end(); tempit++)
	{
		if (tempit->second != a.QDim.at(tempit->first))
		{
			std::cout << "the Q dimention doesn't match each other!" << std::endl;
			exit(1);
		}

	}



	//test the RLQ.
	

	sum.QDim = QDim;
	sum.RLQ = RLQ;
	sum.QMat = QMat;




	for (auto tempit = a.QMat.begin(); tempit != a.QMat.end(); tempit++)
	{
		auto it = QMat.find(tempit->first);
		if (it == QMat.end())
		{
			sum.QMat.insert(std::pair<int, MatrixXd>(tempit->first, a.QMat.at(tempit->first)));
			sum.RLQ.insert(std::pair<int, int>(tempit->first, a.RLQ.at(tempit->first)));
		}

		else
		{
			if (RLQ.at(tempit->first) != a.RLQ.at(tempit->first))
			{
				std::cout << "+: RQ != LQ, WRONG!" << std::endl;
				exit(1);
			}
			sum.QMat.at(tempit->first) += a.QMat.at(tempit->first);
		}
	}




	return sum;
}




void OP::time(const double& x)
{
	for (auto it = QMat.begin(); it != QMat.end(); it++)
	{
		QMat.at(it->first) = QMat.at(it->first) * x;
	}
}







void OP::time(const OP& a, const double& x)
{
	clear();
	QDim = a.QDim;
	RLQ = a.RLQ;
	for (auto it = a.QMat.begin(); it != a.QMat.end(); it++)
	{
		QMat[it->first] = a.QMat.at(it->first)*x;
	}
}



void OP::time(const double& x, const OP& a)
{
	clear();
	QDim = a.QDim;
	RLQ = a.RLQ;
	for (auto it = a.QMat.begin(); it != a.QMat.end(); it++)
	{
		QMat[it->first] = a.QMat.at(it->first)*x;
	}
}






OP OP::operator*(const double& x)
{
	OP times;


	times.QDim = QDim;
	times.RLQ = RLQ;


	//the matrix times the parameter number.
	for (auto tempit = QMat.begin(); tempit != QMat.end(); tempit++)
	{
		times.QMat.insert(std::pair<int, MatrixXd>(tempit->first, x * QMat.at(tempit->first)));
	}

	return times;
}

OP operator*(const double& x, OP& a)
{
	return a*x;
}





void OP::time(const OP&a, const OP&b)
{
	if (a.QDim.size() != b.QDim.size())
	{
		std::cout << "the Q dim size doesn't match!" << std::endl;
		exit(1);
	}
	for (auto tempit = a.QDim.begin(); tempit != a.QDim.end(); tempit++)
	{
		if (tempit->second != b.QDim.at(tempit->first))
		{
			std::cout << "the Q dimention doesn't match each other!" << std::endl;
			exit(1);
		}

	}

	QDim = a.QDim;




	//for the times of two operator, it is a good block systerm. just cross it.
	for (auto ita = b.QMat.begin(); ita != b.QMat.end(); ita++)
	{
		auto it = a.QMat.find(b.RLQ.at(ita->first));
		if (it != a.QMat.end())
		{

			RLQ.insert(std::pair<int, int>(ita->first, a.RLQ.at(it->first)));
			QMat.insert(std::pair<int, MatrixXd>(ita->first, it->second * ita->second));

		}
	}

}




//ltime
int OP::ltime(const OP& a, const OP& wave)
{
	clear();
	int flag(0);
	for (auto ita = wave.QMat.begin(); ita != wave.QMat.end(); ita++)
	{
		auto it = a.QMat.find(wave.RLQ.at(ita->first));
		if (it != a.QMat.end())
		{
			flag++;
			RLQ.insert(std::pair<int, int>(ita->first, a.RLQ.at(it->first)));
			QMat.insert(std::pair<int, MatrixXd>(ita->first, it->second * ita->second));

		}
	}
	return flag;

}






//rtime
int OP::rtime(const OP& a, const OP& wave)
{
	clear();
	int flag(0);

	for (auto ita = a.QMat.begin(); ita != a.QMat.end(); ita++)
	{
		for (auto itQ = wave.RLQ.begin(); itQ != wave.RLQ.end(); itQ++)
		{

			if (ita->first == itQ->first)
			{
				flag++;
				int tempQ(a.RLQ.at(ita->first));
				RLQ[tempQ] = itQ->second;
				QMat[tempQ] = (wave.QMat.at(itQ->first)) * ita->second.transpose();

			}
		}
	}
	return flag;
}





//wave operator, for which don't have the QDim. (the operator in the QWave can't have a QDim)
void OP::addWave(const OP& a, const OP& b)
{
	clear();
	QMat = a.QMat;
	RLQ = a.RLQ;
	for (auto it = b.RLQ.begin(); it != b.RLQ.end(); it++)
	{
		auto tempit = a.RLQ.find(it->first);
		if (tempit != a.RLQ.end())
		{
			if (it->second == tempit->second)
			{
				QMat.at(tempit->first) += b.QMat.at(it->first);
			}
			else
			{
				std::cout << "addWave: a.RLQ.second != b.RLQ.second" << std::endl;
				exit(1);
			}
		}
		else
		{
			RLQ.insert(std::pair<int, int>(it->first, it->second));
			QMat.insert(std::pair<int, MatrixXd>(it->first, b.QMat.at(it->first)));
		}
	}
}




void OP::addWave(const OP& a)
{
	for (auto it = a.RLQ.begin(); it != a.RLQ.end(); it++)
	{
		auto tempit = RLQ.find(it->first);

		if (tempit != RLQ.end())
		{
			QMat.at(tempit->first) += a.QMat.at(it->first);
		}
		else
		{
			RLQ.insert(std::pair<int, int>(it->first, it->second));
			QMat.insert(std::pair<int, MatrixXd>(it->first, a.QMat.at(it->first)));
		}
	}
}






//================for the truncate=======================
void OP::getTruncU(const Parameter& para, const OP& OPWave)
{
	clear();
	std::vector<Eigstruct> Denmat;

	//get the denmat
	for (auto it = OPWave.QMat.begin(); it != OPWave.QMat.end(); it++)
	{
		JacobiSVD<MatrixXd> svd(it->second, ComputeThinU | ComputeThinV);
		for (int i = 0; i<svd.singularValues().size(); i++)
		{
			Eigstruct temp = { OPWave.RLQ.at(it->first), svd.singularValues()(i), svd.matrixU().col(i) };
			Denmat.push_back(temp);
		}
	}
	//std::cout<<"Denmat.size="<<Denmat.size()<<std::endl;
	sort(Denmat.begin(), Denmat.end(), comp);


	int min = Denmat.size() < para.D ? Denmat.size() : para.D;
	//std::cout<<"min = "<<min<<std::endl;
	//get the RLQ/QDim
	for (int i = 0; i<min; i++)
	{
		auto itt = QDim.find(Denmat.at(i).q);
		if (itt != QDim.end())
		{
			itt->second += 1;
		}
		else
		{
			QDim.insert(std::pair<int, int>(Denmat.at(i).q, 1));
			RLQ.insert(std::pair<int, int>(Denmat.at(i).q, Denmat.at(i).q));
		}
	}

	//get the QMat
	for (auto it = QDim.begin(); it != QDim.end(); it++)
	{
		for (int i = 0; i<min; i++)
		{
			if (Denmat.at(i).q == it->first)
			{
				int L = Denmat.at(i).state.size();
				int R = it->second;
				MatrixXd tempmat(L, R);
				QMat.insert(std::pair<int, MatrixXd>(it->first, tempmat));
				break;
			}
		}
	}

	for (auto it = QDim.begin(); it != QDim.end(); it++)
	{
		int ord(0);
		for (int i = 0; i<min; i++)
		{
			if (Denmat.at(i).q == it->first)
			{
				QMat.at(it->first).col(ord) = Denmat.at(i).state;
				ord++;
			}
		}
	}
	
}







void OP::getTruncUR(const Parameter& para, const OP& OPWave)
{
	clear();
	std::vector<Eigstruct> Denmat;

	//get the denmat
	for (auto it = OPWave.QMat.begin(); it != OPWave.QMat.end(); it++)
	{
		JacobiSVD<MatrixXd> svd(it->second, ComputeThinU | ComputeThinV);
		for (int i = 0; i<svd.singularValues().size(); i++)
		{
			Eigstruct temp = { it->first, svd.singularValues()(i), svd.matrixV().col(i) };
			Denmat.push_back(temp);
		}
	}
	//std::cout<<"Denmat.size="<<Denmat.size()<<std::endl;
	sort(Denmat.begin(), Denmat.end(), comp);
	//if(Denmat.size() < para.D) return 0;
	int min = (Denmat.size() < para.D) ? Denmat.size() : para.D;
	//std::cout<<"min = "<<min<<std::endl;
	//get the RLQ/QDim
	for (int i = 0; i<min; i++)
	{
		auto itt = QDim.find(Denmat.at(i).q);
		if (itt != QDim.end())
		{
			itt->second += 1;
		}
		else
		{
			QDim.insert(std::pair<int, int>(Denmat.at(i).q, 1));
			RLQ.insert(std::pair<int, int>(Denmat.at(i).q, Denmat.at(i).q));
		}
	}

	//get the QMat
	for (auto it = QDim.begin(); it != QDim.end(); it++)
	{
		for (int i = 0; i<min; i++)
		{
			if (Denmat.at(i).q == it->first)
			{
				int L = Denmat.at(i).state.size();
				int R = it->second;
				MatrixXd tempmat(L, R);
				QMat.insert(std::pair<int, MatrixXd>(it->first, tempmat));
				break;
			}
		}
	}

	for (auto it = QDim.begin(); it != QDim.end(); it++)
	{
		int ord(0);
		for (int i = 0; i<min; i++)
		{
			if (Denmat.at(i).q == it->first)
			{
				QMat.at(it->first).col(ord) = Denmat.at(i).state;
				ord++;
			}
		}
	}

	//return 1;
}






void OP::getTruncUR(const Parameter& para, const OP& OPWave, double& trance, double& truncerr)
{
	clear();
	std::vector<Eigstruct> Denmat;

	//get the denmat
	for (auto it = OPWave.QMat.begin(); it != OPWave.QMat.end(); it++)
	{
		JacobiSVD<MatrixXd> svd(it->second, ComputeThinU | ComputeThinV);
		for (int i = 0; i<svd.singularValues().size(); i++)
		{
			Eigstruct temp = { it->first, svd.singularValues()(i), svd.matrixV().col(i) };
			Denmat.push_back(temp);
		}
	}
	//std::cout<<"Denmat.size="<<Denmat.size()<<std::endl;
	sort(Denmat.begin(), Denmat.end(), comp);
	//if(Denmat.size() < para.D) return 0;
	int min = (Denmat.size() < para.D) ? Denmat.size() : para.D;
	//std::cout<<"min = "<<min<<std::endl;
	//get the RLQ/QDim
	for (int i = 0; i<min; i++)
	{
		auto itt = QDim.find(Denmat.at(i).q);
		if (itt != QDim.end())
		{
			itt->second += 1;
		}
		else
		{
			QDim.insert(std::pair<int, int>(Denmat.at(i).q, 1));
			RLQ.insert(std::pair<int, int>(Denmat.at(i).q, Denmat.at(i).q));
		}
	}

	//get the QMat
	for (auto it = QDim.begin(); it != QDim.end(); it++)
	{
		for (int i = 0; i<min; i++)
		{
			if (Denmat.at(i).q == it->first)
			{
				int L = Denmat.at(i).state.size();
				int R = it->second;
				MatrixXd tempmat(L, R);
				QMat.insert(std::pair<int, MatrixXd>(it->first, tempmat));
				break;
			}
		}
	}

	for (auto it = QDim.begin(); it != QDim.end(); it++)
	{
		int ord(0);
		for (int i = 0; i<min; i++)
		{
			if (Denmat.at(i).q == it->first)
			{
				QMat.at(it->first).col(ord) = Denmat.at(i).state;
				ord++;
			}
		}
	}

	//return 1;


	trance = 0;
	truncerr = 0;
	for (int i = 0; i<Denmat.size(); i++)
	{
		trance += (Denmat.at(i).lamda)*(Denmat.at(i).lamda);
	}
	for (int i = 0; i< min; i++)
	{
		truncerr += Denmat.at(i).lamda*(Denmat.at(i).lamda);
	}

	truncerr = trance - truncerr;
}







void OP::getTruncU(const Parameter& para, const OP& OPWave, double& trance, double& truncerr)
{
	clear();
	std::vector<Eigstruct> Denmat;

	//get the denmat
	for (auto it = OPWave.QMat.begin(); it != OPWave.QMat.end(); it++)
	{
		JacobiSVD<MatrixXd> svd(it->second, ComputeThinU | ComputeThinV);
		for (int i = 0; i<svd.singularValues().size(); i++)
		{
			Eigstruct temp = { OPWave.RLQ.at(it->first), svd.singularValues()(i), svd.matrixU().col(i) };
			Denmat.push_back(temp);
		}
	}
	//std::cout<<"Denmat.size="<<Denmat.size()<<std::endl;
	sort(Denmat.begin(), Denmat.end(), comp);


	int min = Denmat.size() < para.D ? Denmat.size() : para.D;
	//std::cout<<"min = "<<min<<std::endl;
	//get the RLQ/QDim
	for (int i = 0; i<min; i++)
	{
		auto itt = QDim.find(Denmat.at(i).q);
		if (itt != QDim.end())
		{
			itt->second += 1;
		}
		else
		{
			QDim.insert(std::pair<int, int>(Denmat.at(i).q, 1));
			RLQ.insert(std::pair<int, int>(Denmat.at(i).q, Denmat.at(i).q));
		}
	}

	//get the QMat
	for (auto it = QDim.begin(); it != QDim.end(); it++)
	{
		for (int i = 0; i<min; i++)
		{
			if (Denmat.at(i).q == it->first)
			{
				int L = Denmat.at(i).state.size();
				int R = it->second;
				MatrixXd tempmat(L, R);
				QMat.insert(std::pair<int, MatrixXd>(it->first, tempmat));
				break;
			}
		}
	}

	for (auto it = QDim.begin(); it != QDim.end(); it++)
	{
		int ord(0);
		for (int i = 0; i<min; i++)
		{
			if (Denmat.at(i).q == it->first)
			{
				QMat.at(it->first).col(ord) = Denmat.at(i).state;
				ord++;
			}
		}
	}


	trance = 0;
	truncerr = 0;
	for (int i = 0; i<Denmat.size(); i++)
	{
		trance += pow(Denmat.at(i).lamda, 2);
	}
	for (int i = 0; i< min; i++)
	{
		truncerr += pow(Denmat.at(i).lamda, 2);
	}

	truncerr = trance - truncerr;

}





void OP::DengetTruncU(const Parameter& para, const OP& OPWave, double& trance, double& truncerr)
{
	clear();
	std::vector<Eigstruct> Denmat;

	//get the denmat
	for (auto it = OPWave.QMat.begin(); it != OPWave.QMat.end(); it++)
	{
		SelfAdjointEigenSolver<MatrixXd> es(it->second);
		for (int i = 0; i<es.eigenvalues().size(); i++)
		{
			Eigstruct temp = { it->first, es.eigenvalues()(i), es.eigenvectors().col(i) };
			Denmat.push_back(temp);
		}
	}

	std::stable_sort(Denmat.begin(), Denmat.end(), comp);

	int min = (Denmat.size() < para.D) ? Denmat.size() : para.D;

	//get the RLQ/QDim
	for (int i = 0; i<min; i++)
	{
		auto itt = QDim.find(Denmat.at(i).q);
		if (itt != QDim.end())
		{
			itt->second += 1;
		}
		else
		{
			QDim.insert(std::pair<int, int>(Denmat.at(i).q, 1));
			RLQ.insert(std::pair<int, int>(Denmat.at(i).q, Denmat.at(i).q));
		}
	}

	//get the QMat
	for (auto it = QDim.begin(); it != QDim.end(); it++)
	{
		for (int i = 0; i<min; i++)
		{
			if (Denmat.at(i).q == it->first)
			{
				int L = Denmat.at(i).state.size();
				int R = it->second;
				MatrixXd tempmat(L, R);
				QMat.insert(std::pair<int, MatrixXd>(it->first, tempmat));
				break;
			}
		}
	}

	for (auto it = QDim.begin(); it != QDim.end(); it++)
	{
		int ord(0);
		for (int i = 0; i<min; i++)
		{
			if (Denmat.at(i).q == it->first)
			{
				QMat.at(it->first).col(ord) = Denmat.at(i).state;
				ord++;
			}
		}
	}
	trance = 0;
	truncerr = 0;
	for (int i = 0; i<Denmat.size(); i++)
	{
		trance += Denmat.at(i).lamda;
	}
	for (int i = 0; i< min; i++)
	{
		truncerr += Denmat.at(i).lamda;
	}

	truncerr = trance - truncerr;
}





void OP::DengetTruncU(const Parameter& para, const OP& OPWave, double& trance, double& truncerr, double& Entanglement)
{
        clear();
        std::vector<Eigstruct> Denmat;

        //get the denmat
        for (auto it = OPWave.QMat.begin(); it != OPWave.QMat.end(); it++)
        {
                SelfAdjointEigenSolver<MatrixXd> es(it->second);
                for (int i = 0; i<es.eigenvalues().size(); i++)
                {
                        Eigstruct temp = { it->first, es.eigenvalues()(i), es.eigenvectors().col(i) };
                        Denmat.push_back(temp);
                }
        }

        std::stable_sort(Denmat.begin(), Denmat.end(), comp);

        int min = (Denmat.size() < para.D) ? Denmat.size() : para.D;

        //get the RLQ/QDim
        for (int i = 0; i<min; i++)
        {
                auto itt = QDim.find(Denmat.at(i).q);
                if (itt != QDim.end())
                {
                        itt->second += 1;
                }
                else
                {
                        QDim.insert(std::pair<int, int>(Denmat.at(i).q, 1));
                        RLQ.insert(std::pair<int, int>(Denmat.at(i).q, Denmat.at(i).q));
                }
        }

        //get the QMat
        for (auto it = QDim.begin(); it != QDim.end(); it++)
        {
                for (int i = 0; i<min; i++)
                {
                        if (Denmat.at(i).q == it->first)
                        {
                                int L = Denmat.at(i).state.size();
                                int R = it->second;
                                MatrixXd tempmat(L, R);
                                QMat.insert(std::pair<int, MatrixXd>(it->first, tempmat));
                                break;
                        }
                }
        }

        for (auto it = QDim.begin(); it != QDim.end(); it++)
        {
                int ord(0);
                for (int i = 0; i<min; i++)
                {
                        if (Denmat.at(i).q == it->first)
                        {
                                QMat.at(it->first).col(ord) = Denmat.at(i).state;
                                ord++;
                        }
                }
        }
        trance = 0;
        truncerr = 0;
        Entanglement = 0;
        for (int i = 0; i<Denmat.size(); i++)
        {
                trance += Denmat.at(i).lamda;
                if(Denmat.at(i).lamda > 0.0000000001)
                Entanglement = Entanglement - Denmat.at(i).lamda*log(Denmat.at(i).lamda)/log(2);
                
        }
        for (int i = 0; i< min; i++)
        {
                truncerr += Denmat.at(i).lamda;
        }

        truncerr = trance - truncerr;
}





void OP::truncL(const OP& trunc, const OP& O)
{
	clear();
	if (O.QDim.size() == 0) exit(1);

	for (auto ita = O.QDim.begin(); ita != O.QDim.end(); ita++)
	{
		auto it = trunc.QDim.find(ita->first);
		if (it != trunc.QDim.end())
		{
			QDim.insert(std::pair<int, int>(ita->first, trunc.QDim.at(ita->first)));
		}
	}


	for (auto ita = O.QMat.begin(); ita != O.QMat.end(); ita++)
	{
		auto it = trunc.QMat.find(O.RLQ.at(ita->first));
		if (it != trunc.QMat.end())
		{

			RLQ.insert(std::pair<int, int>(ita->first, it->first));
			QMat.insert(std::pair<int, MatrixXd>(ita->first, it->second.transpose() * ita->second));
		}
	}
}





void OP::truncR(const OP& O, const OP& trunc)
{
	clear();

	for (auto ita = O.QDim.begin(); ita != O.QDim.end(); ita++)
	{

		auto it = trunc.QDim.find(ita->first);

		if (it != trunc.QDim.end())
		{

			QDim.insert(std::pair<int, int>(ita->first, trunc.QDim.at(ita->first)));

		}
	}

	for (auto ita = O.QMat.begin(); ita != O.QMat.end(); ita++)
	{

		auto it = trunc.QMat.find(ita->first);

		if (it != trunc.QMat.end())
		{

			RLQ.insert(std::pair<int, int>(ita->first, O.RLQ.at(ita->first)));
			
			QMat.insert(std::pair<int, MatrixXd>(ita->first, ita->second * it->second));
		}
	}
}





void OP::trunc(const OP& truncU)
{
	OP temp;

	temp.truncR(*this, truncU);

	truncL(truncU, temp);
}






void OP::getDenS(const OP& OPWave)
{
	clear();
	for (auto it = OPWave.RLQ.begin(); it != OPWave.RLQ.end(); it++)
		RLQ.insert(std::pair<int, int>(it->second, it->second));
	for (auto it = OPWave.QMat.begin(); it != OPWave.QMat.end(); it++)
	{
		MatrixXd temp = it->second * (it->second).transpose();
		QMat.insert(std::pair<int, MatrixXd>(OPWave.RLQ.at(it->first), temp));
	}
}







void OP::getDenE(const OP& OPWave)
{
	clear();
	for (auto it = OPWave.RLQ.begin(); it != OPWave.RLQ.end(); it++)
		RLQ.insert(std::pair<int, int>(it->first, it->first));
	for (auto it = OPWave.QMat.begin(); it != OPWave.QMat.end(); it++)
	{
		MatrixXd temp = (it->second).transpose() * it->second;
		QMat.insert(std::pair<int, MatrixXd>(it->first, temp));
	}
}






//to calculate the *this times operator a.
OP OP::operator *(const OP& a)
{
	OP time;



	//to test the QDim.
	if (QDim.size() != a.QDim.size())
	{
		std::cout << "the Q dim size doesn't match!" << std::endl;
		exit(1);
	}
	for (auto tempit = QDim.begin(); tempit != QDim.end(); tempit++)
	{
		if (tempit->second != a.QDim.at(tempit->first))
		{
			std::cout << "the Q dimention doesn't match each other!" << std::endl;
			exit(1);
		}

	}

	time.QDim = QDim;
	MatrixXd tempmat;




	//for the times of two operator, it is a good block systerm. just cross it.
	for (auto ita = a.QMat.begin(); ita != a.QMat.end(); ita++)
	{
		auto it = QMat.find(a.RLQ.at(ita->first));
		if (it != QMat.end())
		{
			tempmat = it->second * ita->second;
			time.RLQ[ita->first] = RLQ[it->first];
			time.QMat[ita->first] = tempmat;

		}
	}

	return time;

}





void OP::operator=(const OP& a)
{
	clear();
	QDim = a.QDim;
	RLQ = a.RLQ;
	QMat = a.QMat;
}





void OP::clear()
{
	QDim.clear();
	QMat.clear();
	RLQ.clear();
}

void OP::show() const
{
	std::cout << "the QDim: " << std::endl;
	for (auto it = QDim.begin(); it != QDim.end(); it++)
	{
		std::cout << it->first << " => " << it->second << std::endl;
	}
	std::cout << "the RLQ: " << std::endl;
	for (auto it = RLQ.begin(); it != RLQ.end(); it++)
	{
		std::cout << it->first << " => " << it->second << std::endl;
	}
	std::cout << "the QMat: " << std::endl;
	for (auto it = QMat.begin(); it != QMat.end(); it++)
	{
		std::cout << it->first << " => " << it->second.rows()<<"x"<<it->second.cols() << std::endl;
	}
}







//save the operator data.
void OP::save(std::ofstream& outfile)
{
	//std::ofstream outfile(str);
	//save the QDim.

	outfile.precision(30);


	outfile << QDim.size() << std::endl;
	for (auto it = QDim.begin(); it != QDim.end(); it++)
	{
		outfile << it->first << "       " << it->second << "        ";
	}


	//save the RLQ.
	outfile << RLQ.size() << std::endl;
	for (auto it = RLQ.begin(); it != RLQ.end(); it++)
	{
		outfile << it->first << "          " << it->second << "             ";
	}
	outfile << std::endl;



	//save the QMat.
	for (auto it = QMat.begin(); it != QMat.end(); it++)
	{
                outfile << it->first << std::endl;

		outfile << it->second << std::endl;
		/*for(int i=0; i<it->second.n_rows; i++)
		{
		for(int j=0; j<it->second.n_cols; j++)
		{
		outfile <<std::setprecision(20)<< it->second(i,j)<<std::endl;
		}
		}*/

	}


	//outfile.close();
}




//read operator data from the file.
void OP::read(std::ifstream& infile)
{
	//std::ifstream infile(str);
	/*if(!infile.is_open())
	{
	std::cout<<"the file "<< str << " is not exit " << std::endl;
	exit(1);
	}*/

	//read in the QDim.
	int size1;
	infile >> size1;

	int tempQ;
	int tempint;
	for (int i = 0; i < size1; i++)
	{

		infile >> tempQ >> tempint;
		QDim[tempQ] = tempint;
	}

	//read in the RLQ.
	int size2;
	infile >> size2;
	int tempQ1;
	for (int i = 0; i < size2; i++)
	{
		infile >> tempQ >> tempQ1;
		RLQ[tempQ] = tempQ1;
	}




	//read in the QMat.
	for (int it = 0; it < RLQ.size(); ++it)
	{

                int tempQR;

                infile >> tempQR;

		MatrixXd A(QDim.at(RLQ.at(tempQR)), QDim.at(tempQR));
		for (int i = 0; i < QDim.at(RLQ.at(tempQR)); i++)
		{
			for (int j = 0; j < QDim.at(tempQR); j++)
			{
				infile >> A(i, j);

			}
		}

		QMat[tempQR] = A;
	}
	//infile.close();



}








//save the truncation operator data.
void OP::truncsave(const int& orbital)
{
        //std::ofstream outfile(str);
        //save the QDim.
        std::string str = itos(orbital);

        str = "./data/trunc" +str;


        std::ofstream outfile(str);

        outfile.precision(30);


        outfile << QDim.size() << std::endl;
        for (auto it = QDim.begin(); it != QDim.end(); it++)
        {
                outfile << it->first << "       " << it->second << "        "<<std::endl;
        }


        //save the RLQ.
        outfile << RLQ.size() << std::endl;
        for (auto it = RLQ.begin(); it != RLQ.end(); it++)
        {
                outfile << it->first << "          " << it->second << "             "<<std::endl;
        }
        outfile << std::endl;



        //save the QMat.
        for (auto it = QMat.begin(); it != QMat.end(); it++)
        {
                outfile << it->first << std::endl;

                outfile << it->second.rows()<<std::endl;

                outfile << it->second.cols() << std::endl;

                outfile << it->second << std::endl;
                

        }


        outfile.close();
}




//read the truncation operator data from the file.
void OP::truncread(const int& orbital)
{
        std::string str = itos(orbital);

        str = "./data/trunc" +str;

        std::ifstream infile(str);
        int size1;
        infile >> size1;

        int tempQ;
        int tempint;
        for (int i = 0; i < size1; i++)
        {

                infile >> tempQ >> tempint;
                QDim[tempQ] = tempint;
        }

        //read in the RLQ.
        int size2;
        infile >> size2;
        int tempQ1;
        for (int i = 0; i < size2; i++)
        {
                infile >> tempQ >> tempQ1;
                RLQ[tempQ] = tempQ1;
        }




        //read in the QMat.
        for (int it = 0; it < RLQ.size(); ++it)
        {

                int tempQR;

                infile >> tempQR;
                int dimL, dimR;
                infile >>dimL >>dimR;

                MatrixXd A(dimL, dimR);
                for (int i = 0; i < dimL; i++)
                {
                        for (int j = 0; j < dimR; j++)
                        {
                                infile >> A(i, j);

                        }
                }

                QMat[tempQR] = A;
        }
        infile.close();



}