#include "QWave.h"


QWave::QWave(){}



QWave::~QWave(){}




QWave::QWave(const QWave& wave)
{
	WavePart = wave.WavePart;
}



void QWave::initial(const QWave& wave)
{
	for(auto it = wave.WavePart.begin(); it != wave.WavePart.end(); ++it)
	{
		for(auto QMatit = it->second.QMat.begin(); QMatit != it->second.QMat.end(); ++QMatit)
		{
			WavePart.at(std::pair<int, int>(it->first.first, it->first.second)).QMat.at(QMatit->first) = QMatit->second;

		}
	}
}






QWave::QWave(const OP& Sys, const OP& m, const OP& n, const OP& Env, int QTot)
{
	for (auto itm = m.QDim.begin(); itm != m.QDim.end(); itm++)
	{
		for (auto itn = n.QDim.begin(); itn != n.QDim.end(); itn++)
		{
			OP tempOP;
			for (auto its = Sys.QDim.begin(); its != Sys.QDim.end(); its++)
			{

				for (auto ite = Env.QDim.begin(); ite != Env.QDim.end(); ite++)
				{
					int q1, q2, qtot_;
					q1=(itm->first+ itn->first);
					q2=(its->first+ ite->first);
					qtot_=(q1+ q2);
					if (qtot_ == QTot)
					{
						tempOP.RLQ[ite->first] = its->first;
						int dimL = its->second;
						int dimR = ite->second;
						MatrixXd tempM(MatrixXd::Zero(dimL, dimR));
						
						tempOP.QMat[ite->first] = tempM;

					}
				}
			}
			if (tempOP.QMat.size()>0)
			{
				WavePart[std::pair<int, int>(itm->first, itn->first)] = tempOP;
			}

		}
	}
}






void QWave::Initial(const OP& Sys, const OP& m, const OP& n, const OP& Env, int QTot)
{
	for (auto itm = m.QDim.begin(); itm != m.QDim.end(); itm++)
	{
		for (auto itn = n.QDim.begin(); itn != n.QDim.end(); itn++)
		{
			OP tempOP;
			for (auto its = Sys.QDim.begin(); its != Sys.QDim.end(); its++)
			{

				for (auto ite = Env.QDim.begin(); ite != Env.QDim.end(); ite++)
				{
					int q1, q2, qtot_;
					q1=(itm->first + itn->first);
					q2=(its->first + ite->first);
					qtot_=(q1 + q2);
					if (qtot_ == QTot)
					{
						tempOP.RLQ[ite->first] = its->first;
						int dimL = its->second;
						int dimR = ite->second;
						MatrixXd tempM(MatrixXd::Zero(dimL, dimR));
						tempOP.QMat[ite->first] = tempM;

					}
				}
			}
			if (tempOP.QMat.size()>0)
			{
				WavePart[std::pair<int, int>(itm->first, itn->first)] = tempOP;
			}

		}
	}
}







void QWave::setZero()
{
	for (auto it = WavePart.begin(); it != WavePart.end(); it++)
	{
		for (auto itt = (it->second).QMat.begin(); itt != (it->second).QMat.end(); itt++)
		{
			int diml = (itt->second).rows();
			int dimr = (itt->second).cols();
			itt->second = (MatrixXd::Zero(diml, dimr));
		}
	}
}






void QWave::normalize()
{
	double x(0);
	for (auto it = WavePart.begin(); it != WavePart.end(); it++)
	{
		for (auto tempit = it->second.QMat.begin(); tempit != it->second.QMat.end(); tempit++)
		{
			int diml = tempit->second.rows();
			int dimr = tempit->second.cols();
			for (int i = 0; i<diml; i++)
			{
				for (int j = 0; j<dimr; j++)
				{
					x +=  ((tempit->second)(i, j)) * ((tempit->second)(i, j));
				}
			}
		}
	}
	double y = sqrt(x);


	for (auto it = WavePart.begin(); it != WavePart.end(); it++)
	{
		for (auto tempit = it->second.QMat.begin(); tempit != it->second.QMat.end(); tempit++)
		{
			tempit->second /= y;
		}
	}
}






int QWave::getDim() const
{
	int n(0);

	for (auto it = WavePart.begin(); it != WavePart.end(); it++)
	{
		for (auto tempit = it->second.QMat.begin(); tempit != it->second.QMat.end(); tempit++)
		{
			n += (tempit->second.rows())*(tempit->second.cols());
		}
	}

	return n;
}





void QWave::add(const QWave&a, const QWave& b)
{
	WavePart = a.WavePart;
	for (auto it = b.WavePart.begin(); it != b.WavePart.end(); it++)
	{
		auto tempit = a.WavePart.find(it->first);
		if (tempit == a.WavePart.end())
		{
			WavePart.insert(std::pair<std::pair<int, int>, OP>(it->first, it->second));
		}
		else
		{
			OP tempOP;
			tempOP.addWave(it->second, tempit->second);
			WavePart.at(it->first) = tempOP;
		}
	}
}







void QWave::Wave2f(std::vector<double>& f) const
{
	f.clear();
	for (auto it = WavePart.begin(); it != WavePart.end(); it++)
	{
		for (auto QMatit = it->second.QMat.begin(); QMatit != it->second.QMat.end(); ++QMatit)
		{
			//auto QMatit = it->second.QMat.find(tempit);
			//if(QMatit == it->second.QMat.end()) continue;
			int diml = QMatit->second.rows();
			int dimr = QMatit->second.cols();
			for (int i = 0; i<diml; i++)
			{
				for (int j = 0; j < dimr; j++)
				{
					double tempv = (QMatit->second)(i, j);
					f.push_back(tempv);
				}
			}
		}
	}
}







void QWave::f2Wave(const std::vector<double>& f)
{
	int n(0);
	for (auto it = WavePart.begin(); it != WavePart.end(); it++)
	{
		for (auto QMatit = it->second.QMat.begin(); QMatit != it->second.QMat.end(); ++QMatit)
		{
			//auto QMatit = it->second.QMat.find(tempit);
			//if(QMatit == it->second.QMat.end()) continue;
			int diml = QMatit->second.rows();
			int dimr = QMatit->second.cols();
			for (int i = 0; i<diml; i++)
			{
				for (int j = 0; j < dimr; j++)
				{
					(QMatit->second)(i, j) = f.at(n);
					n++;
				}
			}
		}
	}
}



void QWave::f2Wave(const VectorXd& f)
{
	int n(0);
	for (auto it = WavePart.begin(); it != WavePart.end(); it++)
	{
		for (auto tempit = it->second.QMat.begin(); tempit != it->second.QMat.end(); tempit++)
		{
			int diml = tempit->second.rows();
			int dimr = tempit->second.cols();
			for (int i = 0; i<diml; i++)
			{
				for (int j = 0; j < dimr; j++)
				{
					(tempit->second)(i, j) = f(n);
					n++;
				}
			}
		}
	}
}







//translate the QWave to OP
void QWave::Wave2OP(OP& O, const OP& sys, const OP& m, const OP& n, const OP& env) const
{
	OP tempS, tempE;
	std::unordered_map<std::pair<int, int>, int, classcom> startDimS, startDimE;
	tempS.findDim(sys, m, tempS.QDim, startDimS);
	tempE.findDim(env, n, tempE.QDim, startDimE);
	//tempS.show();

	for (auto it = WavePart.begin(); it != WavePart.end(); it++)
	{
		for (auto itt = it->second.RLQ.begin(); itt != it->second.RLQ.end(); itt++)
		{
			int tempSQ, tempEQ;
			tempSQ = (it->first.first + itt->second);
			tempEQ = (it->first.second + itt->first);

			int dimS = tempS.QDim.at(tempSQ);
			int dimE = tempE.QDim.at(tempEQ);

			MatrixXd tempmat(MatrixXd::Zero(dimS, dimE));
			
			

			int startL = startDimS.at(std::pair<int, int>(itt->second, it->first.first));
			int startR = startDimE.at(std::pair<int, int>(itt->first, it->first.second));

			int size_row = it->second.QMat.at(itt->first).rows();
			int size_col = it->second.QMat.at(itt->first).cols();

			




			auto temppp = O.QMat.find(tempEQ);
			if (temppp != O.QMat.end())
			{
				tempmat = temppp->second;
			}
			
			//tempmat.block(startL, startR, size_row, size_col) = it->second.QMat.at(itt->first);

			for (int i = 0; i < size_row; ++i)
			{
				for (int j = 0; j < size_col; ++j)
				{
					tempmat(startL + i, startR + j) = it->second.QMat.at(itt->first)(i, j);
				}
			}

			O.QMat[tempEQ] = tempmat;


			auto tempp = O.RLQ.find(tempEQ);
			if (tempp == O.RLQ.end())
			{
				O.RLQ[tempEQ] = tempSQ;
			}

		}
	}

}






//translate the QWave to OP
void QWave::Wave2OP(OP& O, const OP& sys, const OP& m, const OP& n, const OP& env, const int& way) const
{
	O.clear();
	OP tempS, tempE;
	std::unordered_map<std::pair<int, int>, int, classcom> startDimS, startDimE;
	switch (way)
	{
	case 1://==============add m on the right edge of system=============
	{
		       tempS.findDim(sys, m, tempS.QDim, startDimS);
		       tempE.findDim(n, env, tempE.QDim, startDimE);
		       break;
	}
	case -1://==========add n on the left edge of system=============
	{
			tempS.findDim(n, sys, tempS.QDim, startDimS);
			tempE.findDim(env, m, tempE.QDim, startDimE);

			break;
	}
	}
	for (auto it = WavePart.begin(); it != WavePart.end(); it++)
	{
		for (auto itt = it->second.RLQ.begin(); itt != it->second.RLQ.end(); itt++)
		{
			int tempSQ, tempEQ;


			int startL, startR;
			switch (way)
			{
			case 1://==============add m on the right edge of system=============
			{

				       tempSQ = (it->first.first + itt->second);
				       tempEQ = (it->first.second + itt->first);

				       startL = startDimS.at(std::pair<int, int>(itt->second, it->first.first));
				       startR = startDimE.at(std::pair<int, int>(it->first.second, itt->first));

				       break;


			}
			case -1://==========add n on the left edge of system=============
			{
					tempSQ = (itt->second + it->first.second);
					tempEQ = (itt->first + it->first.first);

					startL = startDimS.at(std::pair<int, int>(it->first.second, itt->second));
					startR = startDimE.at(std::pair<int, int>(itt->first, it->first.first));


					break;

			}

			}


			int dimS = tempS.QDim.at(tempSQ);
			int dimE = tempE.QDim.at(tempEQ);

			MatrixXd tempmat(MatrixXd::Zero(dimS, dimE));
			
			

			

			int size_row = it->second.QMat.at(itt->first).rows();
			int size_col = it->second.QMat.at(itt->first).cols();

			




			auto temppp = O.QMat.find(tempEQ);
			if (temppp != O.QMat.end())
			{
				tempmat = temppp->second;
			}
			
			tempmat.block(startL, startR, size_row, size_col) = it->second.QMat.at(itt->first);
			O.QMat[tempEQ] = tempmat;


			auto tempp = O.RLQ.find(tempEQ);
			if (tempp == O.RLQ.end())
			{
				O.RLQ[tempEQ] = tempSQ;
			}

		}
	}

}




//the reroad of some operator.
void QWave::operator=(const QWave& wave)
{
	WavePart.clear();
	WavePart = wave.WavePart;

}



QWave QWave::operator+(const QWave& wave)
{
	QWave sumwave;
	sumwave.WavePart = WavePart;
	for (auto it = wave.WavePart.begin(); it != wave.WavePart.end(); it++)
	{
		auto tempit = WavePart.find(it->first);
		if (tempit == WavePart.end())
		{
			sumwave.WavePart.insert(std::pair<std::pair<int, int>, OP>(it->first, it->second));
		}
		else
		{
			OP tempOP;
			tempOP.addWave(it->second, tempit->second);
			sumwave.WavePart.at(it->first) = tempOP;
		}
	}
	return sumwave;

}






void QWave::show() const
{
	for (auto it = WavePart.begin(); it != WavePart.end(); it++)
	{
		std::cout << "< " << it->first.first << ", " << it->first.second << "> :" << std::endl;
		it->second.show();
	}
}





void QWave::clear()
{
	WavePart.clear();

}






//==============operator function on the QWave 1============================
//|WavePart> = O |wave>
void QWave::OPWave2New(const QWave& wave, const OP& O, int flag)
{
	switch (flag)
	{
	case 1:
	{
		      OSWave2New(O, wave);
		      break;
	}
	case 2:
	{
		      OMWave2New(O, wave);
		      break;
	}
	case 3:
	{
		      ONWave2New(O, wave);
		      break;
	}
	case 4:
	{
		      OEWave2New(O, wave);
		      break;
	}
	}
}






void QWave::OSWave2New(const OP& O, const QWave& wave)
{
	if (O.RLQ.size() == 0)
	{
		exit(1);
	}

	for (auto it = wave.WavePart.begin(); it != wave.WavePart.end(); it++)
	{
		OP temp;
		if (temp.ltime(O, it->second) > 0)
		{
			WavePart.insert(std::pair<std::pair<int, int>, OP>(it->first, temp));
		}
	}
}






void QWave::OEWave2New(const OP& O, const QWave& wave)
{
	if (O.RLQ.size() == 0)
	{
		exit(1);
	}

	for (auto it = wave.WavePart.begin(); it != wave.WavePart.end(); it++)
	{
		OP temp;
		if (temp.rtime(O, it->second) > 0)
		{
			WavePart.insert(std::pair<std::pair<int, int>, OP>(it->first, temp));
		}
	}
}







void QWave::OMWave2New(const OP& O, const QWave& wave)
{
	if (O.RLQ.size() == 0)
	{
		exit(1);
	}

	auto it = O.RLQ.begin();
	int dQ = it->second - it->first;

	for (auto tempit = wave.WavePart.begin(); tempit != wave.WavePart.end(); tempit++)
	{
		auto itt = O.QMat.find(tempit->first.first);
		if (itt != O.QMat.end())
		{
			double x = (itt->second)(0, 0);
			int tempQ(tempit->first.first + dQ);
			OP temp;
			temp.time(x, tempit->second);
			WavePart.insert(std::pair<std::pair<int, int>, OP>(std::pair<int, int>(tempQ, tempit->first.second), temp));
		}
	}
}







void QWave::ONWave2New(const OP& O, const QWave& wave)
{
	if (O.RLQ.size() == 0)
	{
		exit(1);
	}

	auto it = O.RLQ.begin();
	int dQ(it->second - it->first);

	for (auto tempit = wave.WavePart.begin(); tempit != wave.WavePart.end(); tempit++)
	{
		auto itt = O.QMat.find(tempit->first.second);
		if (itt != O.QMat.end())
		{
			double x = (itt->second)(0, 0);
			int tempQ(tempit->first.second + dQ);
			OP temp;
			temp.time(x, tempit->second);
			WavePart.insert(std::pair<std::pair<int, int>, OP>(std::pair<int, int>(tempit->first.first, tempQ), temp));
		}
	}
}







//==================operator function on the QWave 2====================
//|storewave> = O|wavepart> + |storewave>         ==============why it has the second term +|storewave>?==================
void QWave::OPWave(QWave& wave, const OP& O, int flag) const
{
	switch (flag)
	{
	case 1:
	{
		      OSWave(O, wave);
		      break;
	}
	case 2:
	{
		      OMWave(O, wave);
		      break;
	}
	case 3:
	{
		      ONWave(O, wave);
		      break;
	}
	case 4:
	{
		      OEWave(O, wave);
		      break;
	}
	}
}









void QWave::OSWave(const OP& O, QWave& storewave) const
{
	if (O.RLQ.size() == 0)
	{
		exit(1);
	}

	for (auto it = WavePart.begin(); it != WavePart.end(); it++)
	{
		auto itt = storewave.WavePart.find(it->first);
		if (itt != storewave.WavePart.end())
		{
			OP temp;
			if (temp.ltime(O, it->second) > 0)
			{
				itt->second.addWave(temp);
			}
		}
		else
		{
			OP temp;
			if (temp.ltime(O, it->second)>0)
			{
				storewave.WavePart.insert(std::pair<std::pair<int, int>, OP>(it->first, temp));
			}
		}
	}
}






void QWave::OEWave(const OP& O, QWave& storewave) const
{
	if (O.RLQ.size() == 0)
	{
		exit(1);
	}

	for (auto it = WavePart.begin(); it != WavePart.end(); it++)
	{
		auto itt = storewave.WavePart.find(it->first);
		if (itt != storewave.WavePart.end())
		{
			OP temp;
			if (temp.rtime(O, it->second) > 0)
			{
				itt->second.addWave(temp);
			}
		}
		else
		{
			OP temp;
			if (temp.rtime(O, it->second)>0)
			{
				storewave.WavePart.insert(std::pair<std::pair<int, int>, OP>(it->first, temp));
			}
		}
	}
}








void QWave::OMWave(const OP& O, QWave& storewave) const
{
	if (O.RLQ.size() == 0)
	{
		exit(1);
	}

	auto it = O.RLQ.begin();
	int dQ(it->second - it->first);

	QWave tempQW;
	for (auto tempit = WavePart.begin(); tempit != WavePart.end(); tempit++)
	{
		auto itt = O.QMat.find(tempit->first.first);
		if (itt != O.QMat.end())
		{
			double x = (itt->second)(0, 0);
			int tempQ(tempit->first.first + dQ);
			OP temp;
			temp.time(x, tempit->second);
			tempQW.WavePart.insert(std::pair<std::pair<int, int>, OP>(std::pair<int, int>(tempQ, tempit->first.second), temp));
		}
	}
	storewave = tempQW + storewave;
}








void QWave::ONWave(const OP& O, QWave& storewave) const
{
	if (O.RLQ.size() == 0)
	{
		exit(1);
	}

	auto it = O.RLQ.begin();
	int dQ(it->second - it->first);

	QWave tempQW;
	for (auto tempit = WavePart.begin(); tempit != WavePart.end(); tempit++)
	{
		auto itt = O.QMat.find(tempit->first.second);
		if (itt != O.QMat.end())
		{
			double x = (itt->second)(0, 0);
			int tempQ(tempit->first.second + dQ);
			OP temp;
			temp.time(x, tempit->second);
			tempQW.WavePart.insert(std::pair<std::pair<int, int>, OP>(std::pair<int, int>(tempit->first.first, tempQ), temp));
		}
	}
	storewave = tempQW + storewave;
}








//=========for the initial wave======================
//=========================for the Sys eat site M==================================================================================
void QWave::onestepSM(const QWave& wave, const OP&sys, const OP&m, const OP&Env, const OP& n, const OP& truncSM, const OP& truncEN)
{
	clear();
	std::unordered_map<std::pair<int, int>, int, classcom> startDimS;
	std::unordered_map<int, int> nothingDim;

	sys.findDim(sys, m, nothingDim, startDimS);
	//wave.show();m.show();n.show();
	/*std::cout<<"the kron product of sys and m is "<<std::endl;
	for(auto hahait = nothingDim.begin(); hahait != nothingDim.end(); ++hahait)
	{
		std::cout<<hahait->first<<" => "<<hahait->second <<std::endl;
	}*/

	OP truncU;truncU.transO(truncSM);
	


	for(int labeln = 0; labeln < n.RLQ.size(); ++labeln)
	{
		OP tempop;



		for(int labelm = 0; labelm < m.RLQ.size(); ++labelm)
		{
//===================================this part is important======================================================================================
			auto itt = wave.WavePart.find(std::pair<int, int>(labelm, labeln));
			if(itt == wave.WavePart.end()) continue;
//============================for it some parts of wavepart could be disappear when caculate the ground state====================================
			for(auto it = wave.WavePart.at(std::pair<int, int>(labelm,labeln)).RLQ.begin(); it != wave.WavePart.at(std::pair<int, int>(labelm,labeln)).RLQ.end(); ++it)
			{
				auto a = tempop.RLQ.find(it->first);
				if(a == tempop.RLQ.end())
					tempop.RLQ.insert(std::pair<int, int>(it->first, it->second+labelm));


				auto b = tempop.QMat.find(it->first);
				if(b == tempop.QMat.end())
				{
					int dimL(nothingDim.at(it->second + labelm));
					int dimR(wave.WavePart.at(std::pair<int, int>(labelm,labeln)).QMat.at(it->first).cols());
					int MatL(wave.WavePart.at(std::pair<int, int>(labelm,labeln)).QMat.at(it->first).rows());

					MatrixXd tempmat(MatrixXd::Zero(dimL, dimR));

					int startL(startDimS.at(std::pair<int, int>(it->second, labelm)));
					/*std::cout<<startL<<std::endl;
					std::cout<<"the new matrix is "<<dimL <<"x"<<dimR<<std::endl;
					std::cout<<"the old matrix is "<<MatL <<"x"<<dimR<<std::endl;*/

					tempmat.block(startL, 0, MatL, dimR) =  wave.WavePart.at(std::pair<int, int>(labelm,labeln)).QMat.at(it->first);
					tempop.QMat.insert(std::pair<int, MatrixXd>(it->first, tempmat));
				}else
				{
					int dimL(nothingDim.at(it->second + labelm));
					int dimR(wave.WavePart.at(std::pair<int, int>(labelm,labeln)).QMat.at(it->first).cols());
					int MatL(wave.WavePart.at(std::pair<int, int>(labelm,labeln)).QMat.at(it->first).rows());


					int startL(startDimS.at(std::pair<int, int>(it->second, labelm)));
					/*std::cout<<startL<<std::endl;
					std::cout<<"the new matrix is "<<dimL <<"x"<<dimR<<std::endl;
					std::cout<<"the old matrix is "<<MatL <<"x"<<dimR<<std::endl;*/
					b->second.block(startL, 0, MatL, dimR) =  wave.WavePart.at(std::pair<int, int>(labelm,labeln)).QMat.at(it->first);
				}
			}
		}

		
		OP Fop, Ffop;
		//truncU.show();
		//tempop.show();
		Fop.ltime(truncU, tempop);
		Ffop.rtime(truncEN, Fop); 
		//Fop.show();
		//Ffop.show();

		WavePart[std::pair<int, int>(0,labeln)] = Ffop;


	}
}


//==============================fot the sys eats site N================================================================
void QWave::onestepSN(const QWave& wave, const OP&sys, const OP&m, const OP&Env, const OP& n, const OP& truncSN, const OP& truncEM)
{
	clear();
	std::unordered_map<std::pair<int, int>, int, classcom> startDimS;
	std::unordered_map<int, int> nothingDim;

	sys.findDim(n, sys, nothingDim, startDimS);

	OP truncU;truncU.transO(truncSN);
	

	/*wave.show();
	std::cout<<"the kron product of sys and m is "<<std::endl;
	for(auto hahait = nothingDim.begin(); hahait != nothingDim.end(); ++hahait)
	{
		std::cout<<hahait->first<<" => "<<hahait->second <<std::endl;
	}*/

	for(int labelm = 0; labelm < m.RLQ.size(); ++labelm)
	{
		OP tempop;



		for(int labeln = 0; labeln < n.RLQ.size(); ++labeln)
		{
//===================================this part is important======================================================================================
			auto itt = wave.WavePart.find(std::pair<int, int>(labelm, labeln));
			if(itt == wave.WavePart.end()) continue;
//============================for it some parts of wavepart could be disappear when caculate the ground state====================================
			for(auto it = wave.WavePart.at(std::pair<int, int>(labelm,labeln)).RLQ.begin(); it != wave.WavePart.at(std::pair<int, int>(labelm,labeln)).RLQ.end(); ++it)
			{
				auto a = tempop.RLQ.find(it->first);
				if(a == tempop.RLQ.end())
					tempop.RLQ.insert(std::pair<int, int>(it->first, it->second+labeln));


				auto b = tempop.QMat.find(it->first);
				if(b == tempop.QMat.end())
				{
					int dimL(nothingDim.at(it->second + labeln));
					int dimR(wave.WavePart.at(std::pair<int, int>(labelm,labeln)).QMat.at(it->first).cols());
					int MatL(wave.WavePart.at(std::pair<int, int>(labelm,labeln)).QMat.at(it->first).rows());

					MatrixXd tempmat(MatrixXd::Zero(dimL, dimR));

					int startL(startDimS.at(std::pair<int, int>(labeln, it->second)));
					/*std::cout<<startL<<std::endl;
					std::cout<<"the new matrix is "<<dimL <<"x"<<dimR<<std::endl;
					std::cout<<"the old matrix is "<<MatL <<"x"<<dimR<<std::endl;*/
					tempmat.block(startL, 0, MatL, dimR) =  wave.WavePart.at(std::pair<int, int>(labelm,labeln)).QMat.at(it->first);
					tempop.QMat.insert(std::pair<int, MatrixXd>(it->first, tempmat));
				}else
				{
					int dimL(nothingDim.at(it->second + labeln));
					int dimR(wave.WavePart.at(std::pair<int, int>(labelm,labeln)).QMat.at(it->first).cols());
					int MatL(wave.WavePart.at(std::pair<int, int>(labelm,labeln)).QMat.at(it->first).rows());

					int startL(startDimS.at(std::pair<int, int>(labeln, it->second)));
					/*std::cout<<startL<<std::endl;
					std::cout<<"the new matrix is "<<dimL <<"x"<<dimR<<std::endl;
					std::cout<<"the old matrix is "<<MatL <<"x"<<dimR<<std::endl;*/
					b->second.block(startL, 0, MatL, dimR) =  wave.WavePart.at(std::pair<int, int>(labelm,labeln)).QMat.at(it->first);
				}
			}
		}

		
		OP Fop, Ffop;
		//truncU.show();
		//tempop.show();
		Fop.ltime(truncU, tempop);
		Ffop.rtime(truncEM, Fop);
		WavePart[std::pair<int, int>(labelm,0)] = Ffop;


	}
}



//==================================for the Env eat the site M===========================================================
void QWave::onestepEM(const QWave& wave, const OP&sys, const OP&m, const OP&Env, const OP& n, const OP& truncEM, const OP& truncSN)
{
	clear();
	std::unordered_map<std::pair<int, int>, int, classcom> startDimS;
	std::unordered_map<int, int> nothingDim;

	sys.findDim(Env, m, nothingDim, startDimS);
	//wave.show();m.show();n.show();
	/*std::cout<<"the kron product of sys and m is "<<std::endl;
	for(auto hahait = nothingDim.begin(); hahait != nothingDim.end(); ++hahait)
	{
		std::cout<<hahait->first<<" => "<<hahait->second <<std::endl;
	}*/

	OP truncV;truncV.transO(truncEM);
	


	for(int labeln = 0; labeln < n.RLQ.size(); ++labeln)
	{
		OP tempop;



		for(int labelm = 0; labelm < m.RLQ.size(); ++labelm)
		{
//===================================this part is important======================================================================================
			auto itt = wave.WavePart.find(std::pair<int, int>(labelm, labeln));
			if(itt == wave.WavePart.end()) continue;
//============================for it some parts of wavepart could be disappear when caculate the ground state====================================
			for(auto it = wave.WavePart.at(std::pair<int, int>(labelm,labeln)).RLQ.begin(); it != wave.WavePart.at(std::pair<int, int>(labelm,labeln)).RLQ.end(); ++it)
			{
				auto a = tempop.RLQ.find(it->first + labelm);
				if(a == tempop.RLQ.end())
					tempop.RLQ.insert(std::pair<int, int>(it->first + labelm, it->second));


				auto b = tempop.QMat.find(it->first + labelm);
				if(b == tempop.QMat.end())
				{
					int dimL(nothingDim.at(it->first + labelm));
					int dimR(wave.WavePart.at(std::pair<int, int>(labelm,labeln)).QMat.at(it->first).cols());
					int MatL(wave.WavePart.at(std::pair<int, int>(labelm,labeln)).QMat.at(it->first).rows());

					MatrixXd tempmat(MatrixXd::Zero(MatL, dimL));

					int startR(startDimS.at(std::pair<int, int>(it->first, labelm)));
					/*std::cout<<startL<<std::endl;
					std::cout<<"the new matrix is "<<dimL <<"x"<<dimR<<std::endl;
					std::cout<<"the old matrix is "<<MatL <<"x"<<dimR<<std::endl;*/

					tempmat.block(0, startR, MatL, dimR) =  wave.WavePart.at(std::pair<int, int>(labelm,labeln)).QMat.at(it->first);
					tempop.QMat.insert(std::pair<int, MatrixXd>(it->first + labelm, tempmat));
				}else
				{
					int dimL(nothingDim.at(it->first + labelm));
					int dimR(wave.WavePart.at(std::pair<int, int>(labelm,labeln)).QMat.at(it->first).cols());
					int MatL(wave.WavePart.at(std::pair<int, int>(labelm,labeln)).QMat.at(it->first).rows());


					int startR(startDimS.at(std::pair<int, int>(it->first, labelm)));
					/*std::cout<<startL<<std::endl;
					std::cout<<"the new matrix is "<<dimL <<"x"<<dimR<<std::endl;
					std::cout<<"the old matrix is "<<MatL <<"x"<<dimR<<std::endl;*/
					b->second.block(0, startR, MatL, dimR) =  wave.WavePart.at(std::pair<int, int>(labelm,labeln)).QMat.at(it->first);
				}
			}
		}

		
		OP Fop, Ffop;
		//truncEM.show();
		//tempop.show();
		Fop.ltime(truncSN, tempop);
		Ffop.rtime(truncV, Fop); 
		WavePart[std::pair<int, int>(0,labeln)] = Ffop;


	}
}


//=============================for the Env eat the site N==========================================================================
void QWave::onestepEN(const QWave& wave, const OP&sys, const OP&m, const OP&Env, const OP& n, const OP& truncEN, const OP& truncSM)
{
	clear();
	std::unordered_map<std::pair<int, int>, int, classcom> startDimS;
	std::unordered_map<int, int> nothingDim;

	sys.findDim(n, Env, nothingDim, startDimS);

	OP truncV;truncV.transO(truncEN);
	

	/*wave.show();
	std::cout<<"the kron product of sys and m is "<<std::endl;
	for(auto hahait = nothingDim.begin(); hahait != nothingDim.end(); ++hahait)
	{
		std::cout<<hahait->first<<" => "<<hahait->second <<std::endl;
	}*/

	for(int labelm = 0; labelm < m.RLQ.size(); ++labelm)
	{
		OP tempop;



		for(int labeln = 0; labeln < n.RLQ.size(); ++labeln)
		{
//===================================this part is important======================================================================================
			auto itt = wave.WavePart.find(std::pair<int, int>(labelm, labeln));
			if(itt == wave.WavePart.end()) continue;
//============================for it some parts of wavepart could be disappear when caculate the ground state====================================
			for(auto it = wave.WavePart.at(std::pair<int, int>(labelm,labeln)).RLQ.begin(); it != wave.WavePart.at(std::pair<int, int>(labelm,labeln)).RLQ.end(); ++it)
			{
				auto a = tempop.RLQ.find(it->first + labeln);
				if(a == tempop.RLQ.end())
					tempop.RLQ.insert(std::pair<int, int>(it->first + labeln, it->second));


				auto b = tempop.QMat.find(it->first + labeln);
				if(b == tempop.QMat.end())
				{
					int dimL(nothingDim.at(it->first + labeln));
					int dimR(wave.WavePart.at(std::pair<int, int>(labelm,labeln)).QMat.at(it->first).cols());
					int MatL(wave.WavePart.at(std::pair<int, int>(labelm,labeln)).QMat.at(it->first).rows());

					MatrixXd tempmat(MatrixXd::Zero(MatL, dimL));

					int startR(startDimS.at(std::pair<int, int>(labeln, it->first)));
					/*std::cout<<startL<<std::endl;
					std::cout<<"the new matrix is "<<dimL <<"x"<<dimR<<std::endl;
					std::cout<<"the old matrix is "<<MatL <<"x"<<dimR<<std::endl;*/
					tempmat.block(0, startR, MatL, dimR) =  wave.WavePart.at(std::pair<int, int>(labelm,labeln)).QMat.at(it->first);
					tempop.QMat.insert(std::pair<int, MatrixXd>(it->first + labeln, tempmat));
				}else
				{
					int dimL(nothingDim.at(it->first + labeln));
					int dimR(wave.WavePart.at(std::pair<int, int>(labelm,labeln)).QMat.at(it->first).cols());
					int MatL(wave.WavePart.at(std::pair<int, int>(labelm,labeln)).QMat.at(it->first).rows());

					int startR(startDimS.at(std::pair<int, int>(labeln, it->first)));
					/*std::cout<<startL<<std::endl;
					std::cout<<"the new matrix is "<<dimL <<"x"<<dimR<<std::endl;
					std::cout<<"the old matrix is "<<MatL <<"x"<<dimR<<std::endl;*/
					b->second.block(0, startR, MatL, dimR) =  wave.WavePart.at(std::pair<int, int>(labelm,labeln)).QMat.at(it->first);
				}
			}
		}

		
		OP Fop, Ffop;
		//truncEN.show();
		//tempop.show();
		Fop.ltime(truncSM, tempop);
		Ffop.rtime(truncV, Fop);
		WavePart[std::pair<int, int>(labelm,0)] = Ffop;


	}
}



//===========================after the Sys eat site M, the Env spit site N====================
void QWave::twostepSM(const QWave&wave, const OP&Sys, const OP&m, const OP& Env, const OP&n)
{
	clear();
	std::unordered_map<std::pair<int, int>, int, classcom> startDim;
	std::unordered_map<int, int> nothingDim;

	Sys.findDim(Env, m, nothingDim, startDim);
	//wave.show();
	for(auto labeln = wave.WavePart.begin(); labeln != wave.WavePart.end(); ++labeln)
	{

		for(int labelm = 0; labelm < m.QDim.size(); ++labelm)
		{
			OP tempop;
			for(auto it = labeln->second.QMat.begin(); it != labeln->second.QMat.end(); ++it)
			{
				int EnvN(it->first - labelm);
				auto tempit = startDim.find(std::pair<int, int>(EnvN, labelm));
				if(tempit == startDim.end()) continue;

				auto RLQit = tempop.RLQ.find(EnvN);
				if(RLQit == tempop.RLQ.end()) 
					tempop.RLQ.insert(std::pair<int, int>(EnvN, labeln->second.RLQ.at(it->first)));

				int DimL(it->second.rows());
				int DimR(Env.QMat.at(EnvN).cols());
				int start(tempit->second);

				auto QMatit = tempop.QMat.find(EnvN);
				if(QMatit == tempop.QMat.end())
				{
					//std::cout<<"the old dim: "<<it->second.rows()<<"x"<<it->second.cols()<<std::endl
					//<<"the new dim: "<<DimL<<"x"<<DimR<<" begin at "<<"("<<0<<","<<start<<")"<<std::endl;

					MatrixXd tempmat(it->second.block(0, start, DimL, DimR));


					tempop.QMat.insert(std::pair<int, MatrixXd>(EnvN, tempmat));

					/*std::cout << "the QMat: "<<std::endl
					<<EnvN<<" => "<<tempop.QMat.at(EnvN);*/
				}

			}
			//std::cout<<"<"<<labelm<<", "<<labeln->first.second<<">"<<std::endl;
			//tempop.show();
			if(tempop.QMat.size() != 0)
			WavePart[std::pair<int, int>(labelm, labeln->first.second)] = tempop;
		}
	}
}





//===========================after the Sys eat site M, the Env spit site N====================
void QWave::twostepSN(const QWave&wave, const OP&Sys, const OP&m, const OP& Env, const OP&n)
{
	clear();
	std::unordered_map<std::pair<int, int>, int, classcom> startDim;
	std::unordered_map<int, int> nothingDim;

	Sys.findDim(n, Env, nothingDim, startDim);
	/*for(auto tempitt = startDim.begin(); tempitt != startDim.end(); ++tempitt)
	{
		std::cout<<"<"<<tempitt->first.first<<", "<<tempitt->first.second<<">:"<<tempitt->second<<std::endl;
	}exit(true);*/
	//wave.show();
	for(auto labelm = wave.WavePart.begin(); labelm != wave.WavePart.end(); ++labelm)
	{

		for(int labeln = 0; labeln < n.QDim.size(); ++labeln)
		{
			OP tempop;
			for(auto it = labelm->second.QMat.begin(); it != labelm->second.QMat.end(); ++it)
			{
				int EnvN(it->first - labeln);
				auto tempit = startDim.find(std::pair<int, int>(labeln, EnvN));
				if(tempit == startDim.end()) continue;

				auto RLQit = tempop.RLQ.find(EnvN);
				if(RLQit == tempop.RLQ.end()) 
					tempop.RLQ.insert(std::pair<int, int>(EnvN, labelm->second.RLQ.at(it->first)));

				int DimL(it->second.rows());
				int DimR(Env.QMat.at(EnvN).cols());
				int start(tempit->second);

				auto QMatit = tempop.QMat.find(EnvN);
				if(QMatit == tempop.QMat.end())
				{
					//std::cout<<"the old dim: "<<it->second.rows()<<"x"<<it->second.cols()<<std::endl
					//<<"the new dim: "<<DimL<<"x"<<DimR<<" begin at "<<"("<<0<<","<<start<<")"<<std::endl;

					MatrixXd tempmat(it->second.block(0, start, DimL, DimR));


					tempop.QMat.insert(std::pair<int, MatrixXd>(EnvN, tempmat));

					/*std::cout << "the QMat: "<<std::endl
					<<EnvN<<" => "<<tempop.QMat.at(EnvN);*/
				}

			}
			//std::cout<<"<"<<labelm<<", "<<labeln->first.second<<">"<<std::endl;
			//tempop.show();
			if(tempop.QMat.size() != 0)
			WavePart[std::pair<int, int>(labelm->first.first, labeln)] = tempop;
		}
	}
}




void QWave::twostepEM(const QWave&wave, const OP&Sys, const OP&m, const OP& Env, const OP&n)
{
	clear();
	std::unordered_map<std::pair<int, int>, int, classcom> startDim;
	std::unordered_map<int, int> nothingDim;

	Sys.findDim(Sys, m, nothingDim, startDim);
	//wave.show();
	for(auto labeln = wave.WavePart.begin(); labeln != wave.WavePart.end(); ++labeln)
	{

		for(int labelm = 0; labelm < m.QDim.size(); ++labelm)
		{
			OP tempop;
			for(auto it = labeln->second.QMat.begin(); it != labeln->second.QMat.end(); ++it)
			{
				int SysN(labeln->second.RLQ.at(it->first) - labelm);
				auto tempit = startDim.find(std::pair<int, int>(SysN, labelm));
				if(tempit == startDim.end()) continue;

				auto RLQit = tempop.RLQ.find(it->first);
				if(RLQit == tempop.RLQ.end()) 
					tempop.RLQ.insert(std::pair<int, int>(it->first, SysN));

				int DimL(Sys.QMat.at(SysN).rows());
				int DimR(it->second.cols());
				int start(tempit->second);

				auto QMatit = tempop.QMat.find(it->first);
				if(QMatit == tempop.QMat.end())
				{
					//std::cout<<"the old dim: "<<it->second.rows()<<"x"<<it->second.cols()<<std::endl
					//<<"the new dim: "<<DimL<<"x"<<DimR<<" begin at "<<"("<<0<<","<<start<<")"<<std::endl;

					MatrixXd tempmat(it->second.block(start, 0, DimL, DimR));


					tempop.QMat.insert(std::pair<int, MatrixXd>(it->first, tempmat));

					/*std::cout << "the QMat: "<<std::endl
					<<EnvN<<" => "<<tempop.QMat.at(EnvN);*/
				}

			}
			//std::cout<<"<"<<labelm<<", "<<labeln->first.second<<">"<<std::endl;
			//tempop.show();
			if(tempop.QMat.size() != 0)
			WavePart[std::pair<int, int>(labelm, labeln->first.second)] = tempop;
		}
	}
}



void QWave::twostepEN(const QWave&wave, const OP&Sys, const OP&m, const OP& Env, const OP&n)
{
	clear();
	std::unordered_map<std::pair<int, int>, int, classcom> startDim;
	std::unordered_map<int, int> nothingDim;

	Sys.findDim(n, Sys, nothingDim, startDim);
	/*for(auto tempitt = startDim.begin(); tempitt != startDim.end(); ++tempitt)
	{
		std::cout<<"<"<<tempitt->first.first<<", "<<tempitt->first.second<<">:"<<tempitt->second<<std::endl;
	}*/
	//wave.show();
	for(auto labelm = wave.WavePart.begin(); labelm != wave.WavePart.end(); ++labelm)
	{

		for(int labeln = 0; labeln < n.QDim.size(); ++labeln)
		{
			OP tempop;
			for(auto it = labelm->second.QMat.begin(); it != labelm->second.QMat.end(); ++it)
			{
				int SysN(labelm->second.RLQ.at(it->first) - labeln);
				auto tempit = startDim.find(std::pair<int, int>(labeln, SysN));
				if(tempit == startDim.end()) continue;

				auto RLQit = tempop.RLQ.find(it->first);
				if(RLQit == tempop.RLQ.end()) 
					tempop.RLQ.insert(std::pair<int, int>(it->first, SysN));

				int DimL(Sys.QMat.at(SysN).rows());
				int DimR(it->second.cols());
				int start(tempit->second);

				auto QMatit = tempop.QMat.find(it->first);
				if(QMatit == tempop.QMat.end())
				{
					/*std::cout<<"the old dim: "<<it->second.rows()<<"x"<<it->second.cols()<<std::endl
					<<"the new dim: "<<DimL<<"x"<<DimR<<" begin at "<<"("<<0<<","<<start<<")"<<std::endl;*/

					MatrixXd tempmat(it->second.block(start, 0, DimL, DimR));


					tempop.QMat.insert(std::pair<int, MatrixXd>(it->first, tempmat));

					/*std::cout << "the QMat: "<<std::endl
					<<EnvN<<" => "<<tempop.QMat.at(EnvN);*/
				}

			}
			//std::cout<<"<"<<labelm<<", "<<labeln->first.second<<">"<<std::endl;
			//tempop.show();
			if(tempop.QMat.size() != 0)
			WavePart[std::pair<int, int>(labelm->first.first, labeln)] = tempop;
		}
	}
}