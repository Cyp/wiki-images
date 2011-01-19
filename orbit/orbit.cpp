#include "../headers/image.h"
#include <eigen2/Eigen/Core>
#include <eigen2/Eigen/LU>
USING_PART_OF_NAMESPACE_EIGEN


static const int D = 2;

struct System
{
	System(int num);
	void init();
	void render(Image &image);
	void updateQ(double dt);
	void updateP(double dt);
	void matMul();
	Pixel &C(int i) { return cVec[i]; }
	double &M(int i) { return mVec[i]; }
	double &Q(int i, int d) { return zVec[(i*D + d)*2 + 0]; }
	double &P(int i, int d) { return zVec[(i*D + d)*2 + 1]; }
	double &Aqq(int i1, int d1, int i2, int d2) { return aMat((i1*D + d1)*2 + 0, (i2*D + d2)*2 + 0); }
	double &Aqp(int i1, int d1, int i2, int d2) { return aMat((i1*D + d1)*2 + 0, (i2*D + d2)*2 + 1); }
	double &Apq(int i1, int d1, int i2, int d2) { return aMat((i1*D + d1)*2 + 1, (i2*D + d2)*2 + 0); }
	double &App(int i1, int d1, int i2, int d2) { return aMat((i1*D + d1)*2 + 1, (i2*D + d2)*2 + 1); }

	int num;
	int params;
	std::vector<Pixel>  cVec;
	std::vector<double> mVec;
	VectorXd            zVec;
	VectorXd            zVecInitial;
	MatrixXd            aMat;
	MatrixXd            aMatTotal;
};

void System::render(Image &image)
{
	for (int i = 0; i < num; ++i)
	{
		double radius = pow(M(i), 1/3.)/100;
		image.addEllipse(Q(i, 0) - radius, Q(i, 1) - radius, Q(i, 0) + radius, Q(i, 1) + radius, -C(i));
	}
}

System::System(int num)
	: num(num)
	, params(num*D*2)
	, zVec(params)
	, aMat(MatrixXd::Identity(params, params))
	, aMatTotal(MatrixXd::Identity(params, params))
{
	for (int i = 0; i < num; ++i)
	{
		cVec.push_back(Pixel(rand()*(1./RAND_MAX), rand()*(1./RAND_MAX), rand()*(1./RAND_MAX)));
	}
	for (int i = 0; i < num; ++i)
	{
		mVec.push_back(0.5); //rand()*(0.1/RAND_MAX));
	}
	for (int i = 0; i < num; ++i)
	{
		for (int d = 0; d < D; ++d)
		{
			Q(i, d) = rand()*(1./RAND_MAX);
			P(i, d) = (-1 + rand()*(2./RAND_MAX)) * 1.5;
		}
	}

	zVecInitial = zVec;
}

void System::init()
{
	zVec = zVecInitial;
	aMat.setIdentity(params, params);
	aMatTotal.setIdentity(params, params);
}

void System::updateQ(double dt)
{
	for (int i = 0; i < num; ++i)
	{
		for (int d = 0; d < D; ++d)
		{
			Q(i, d) += P(i, d) / M(i) * dt;
			Aqp(i, d, i, d) += 1 / M(i) * dt;
		}
	}

	matMul();
}

static inline double sq(double x) { return x*x; }

void System::updateP(double dt)
{
	for (int i = 0; i < num; ++i)
		for (int j = 0; j < i; ++j)
	{
		double distSq = 0;
		for (int e = 0; e < D; ++e)
		{
			distSq += sq(Q(j, e) - Q(i, e)) + 1;
		}
		double f = M(i)*M(j) / (distSq*sqrt(distSq));
		for (int d = 0; d < D; ++d)
		{
			P(i, d) +=  f * (Q(j, d) - Q(i, d)) * dt;
			P(j, d) += -f * (Q(j, d) - Q(i, d)) * dt;
			Apq(i, d, j, d) +=  f * dt;
			Apq(i, d, i, d) += -f * dt;
			Apq(j, d, j, d) += -f * dt;
			Apq(j, d, i, d) +=  f * dt;
			double g = -3/2. * f * (Q(j, d) - Q(i, d)) / distSq;
			for (int e = 0; e < D; ++e)
			{
				double h = 2 * (Q(j, e) - Q(i, e));
				Apq(i, d, j, e) +=  g * h * dt;
				Apq(i, d, i, e) += -g * h * dt;
				Apq(j, d, j, e) += -g * h * dt;
				Apq(j, d, i, e) +=  g * h * dt;
			}
		}
	}

	matMul();
}

void System::matMul()
{
	aMatTotal = aMatTotal*aMat;
	aMat = MatrixXd::Identity(params, params);
}

const int imageSize = 400;
const double border = 2;
void evolve(System &system, bool doRender)
{
	system.init();

	double totalTime = 1;
	int frames = 100;
	double step = 0.01;
	for (int frame = 0; frame < frames; ++frame)
	{
		if (doRender)
		{
			char name[100];
			sprintf(name, "OrbitAnim%04d", frame);
			Image image(imageSize, imageSize, 0 - border, 0 - border, 1 + border, 1 + border);
			system.render(image);
			image.save(name);
		}

		double frameTime = totalTime / frames;
		while (frameTime > 1.e-30)
		{
			double dt = std::min(step, frameTime);
			frameTime -= dt;

			system.updateP(7/24. * dt);
			system.updateQ(2/3. * dt);
			system.updateP(3/4. * dt);
			system.updateQ(-2/3. * dt);
			system.updateP(-1/24. * dt);
			system.updateQ(1. * dt);
		}
	}
}

Eigen::LU<MatrixXd> getLU(MatrixXd const &matrix, MatrixXd const &kernelMatrix)
{
	for (double n = 1; n < 1e100; n *= 1000)
	{
		Eigen::LU<MatrixXd> result = (matrix + n*kernelMatrix).lu();
		if (result.isInvertible())
			return result;
	}
	Eigen::LU<MatrixXd> result = (matrix + 1.e12*kernelMatrix).lu();
	if (!result.isInvertible())
	{
		MatrixXd mat = result.kernel();
		printf("Matrix not invertible!\nMatrix =\n");
		for (int i = 0; i < matrix.rows(); ++i)
		{
			for (int j = 0; j < matrix.cols(); ++j)
			{
				printf("\t%.2lf", matrix(i, j));
			}
			printf("\n");
		}
		printf("Kernel =\n");
		for (int i = 0; i < mat.rows(); ++i)
		{
			for (int j = 0; j < mat.cols(); ++j)
			{
				printf("\t%.2lf", mat(i, j));
			}
			printf("\n");
		}
	}
	return result;
}

VectorXd getAdjustVector(System &system, MatrixXd &matShuffle, MatrixXd &kernelMatrix)
{
	evolve(system, false);
	VectorXd diff = system.zVec - matShuffle*system.zVecInitial;
	VectorXd adjust;
	getLU(system.aMatTotal - matShuffle*system.aMat, kernelMatrix).solve(diff, &adjust);
	//printf("diff = %lf, adjust = %lf\n", diff.norm(), adjust.norm());
	if (adjust.size() == 0)
	{
		return adjust;
		adjust.resize(system.params);
		for (int i = 0; i < system.params; ++i)
			adjust(i) = (-1 + rand()*(2./RAND_MAX));
	}
	//double maxAdjust = 0.1;
	//double adjustAmount = tanh(adjust.norm()/maxAdjust)*maxAdjust;

	return -adjust;// * (adjustAmount/adjust.norm());
}

int main()
{
	srand(13);
	System system(13);
	int shuffleMap[7] = {2, 0, 1, 6, 3, 4, 5};
	system.M(0) = 10;
	/*system.M(0) = 70;
	system.Q(0, 0) = 0;
	system.Q(0, 1) = -1;
	system.P(0, 0) = 100;
	system.P(0, 1) = 0;

	system.M(1) = 70;
	system.Q(1, 0) = 0;
	system.Q(1, 1) = 1;
	system.P(1, 0) = -100;
	system.P(1, 1) = 0;*/
	system.zVecInitial = system.zVec;
	{
		Image image(imageSize, imageSize, 0 - border, 0 - border, 1 + border, 1 + border);
		system.render(image);
		image.save("Orbit");
	}

	MatrixXd matShuffle(system.params, system.params);
	/*for (int i = 0; i < system.params; ++i)
		for (int j = 0; j < system.params; ++j)
			matShuffle(i - i/(D*2)*D*2 + shuffleMap[i/(D*2)]*D*2, j) = i == j;*/

	MatrixXd kernelMatrix(system.params, system.params);
	for (int i = 0; i < system.params; ++i)
		for (int j = 0; j < system.params; ++j)
			kernelMatrix(i, j) = (i%2 == 1 && j%2 == 0 && i/2%D == j/2%D);

	if (true)
	{
		matShuffle = MatrixXd::Identity(system.params, system.params);
	}

	//for (int n = 0; n < 100; ++n)
	//double adjustDiv = 10;
	while (true)
	{
		//adjustDiv -= (adjustDiv - 0.5)*0.01;
		/*char name[100];
		sprintf(name, "OrbitAnim%04d", n);
		Image image(400, 400, 0, 0, 1, 1);
		system.render(image);
		evolve(system, false);
		system.render(image);
		image.save(name);*/

		evolve(system, false);
		//system.updateP(10);

		VectorXd diff = system.zVec - matShuffle*system.zVecInitial;
		if (diff.norm() < 0.1)
		{
			printf("Done!\n");
			system.init();
			break;
		}
/*Eigen::LU<MatrixXd> result = (system.aMatTotal - matShuffle*system.aMat, kernelMatrix).lu();
for (int i = 0; i < system.params; ++i)
{
	printf("\t");
	for (int j = 0; j < system.params; ++j)
	{
		printf("\t%.2lf", (system.aMatTotal - matShuffle*system.aMat, kernelMatrix)(i, j));
	}
	printf("\n");
}
printf("\n");
if (result.isInvertible())
{
MatrixXd mat = result.inverse();
for (int i = 0; i < mat.rows(); ++i)
{
	printf("%.2lf\t", (mat*diff)(i));
	for (int j = 0; j < mat.cols(); ++j)
	{
		printf("\t%.2lf", mat(i, j));
	}
	printf("\t\t%.2lf\n", diff(i));
}
}
else
{
MatrixXd mat = result.kernel();
for (int i = 0; i < mat.rows(); ++i)
{
	printf("?\t");
	for (int j = 0; j < mat.cols(); ++j)
	{
		printf("\t%.2lf", mat(i, j));
	}
	printf("\t\t%.2lf\n", diff(i));
}
}*/


		VectorXd adjust;
		getLU(system.aMatTotal - matShuffle*system.aMat, kernelMatrix).solve(diff, &adjust);
		printf("diff = %lf, adjust = %lf\n", diff.norm(), adjust.norm());
		if (diff.norm() > 10000 || adjust.size() == 0 || system.zVecInitial.norm() > 20)
		{
			printf("Reset!\n");
			system = System(system.num);
			continue;
		}
#if 0
		if (adjust.size() == 0)
		{
			adjust.resize(system.params);
			for (int i = 0; i < system.params; ++i)
				adjust(i) = (-1 + rand()*(2./RAND_MAX));
		}
#endif
#if 0
		double adjustAmount = std::min(1., adjust.norm()/adjustDiv);
		system.zVecInitial -= adjust * (adjustAmount/adjust.norm());
		system.init();
#endif
		double h = 0.01;
		VectorXd tmp = system.zVecInitial;
		system.zVecInitial = tmp;
		VectorXd k1 = getAdjustVector(system, matShuffle, kernelMatrix);
		if (k1.size() == 0) { printf("Reset!\n"); system = System(system.num); continue; }
//for (int i = 0; i < k1.size(); ++i)
//	k1(i) = (-1 + rand()*(2./RAND_MAX));
//k1 *= 1;
//VectorXd oldDiff = system.zVec - matShuffle*system.zVecInitial;
//VectorXd predictedDiffDelta = (system.aMatTotal - matShuffle*system.aMat) * (h * 1/2 * k1);
		system.zVecInitial = tmp + h * 1/2 * k1;
		VectorXd k2 = getAdjustVector(system, matShuffle, kernelMatrix);
		if (k2.size() == 0) { printf("Reset!\n"); system = System(system.num); continue; }
//VectorXd newDiff = system.zVec - matShuffle*system.zVecInitial;
//VectorXd diffDelta = newDiff - diff;
//	printf("Pre ");
//	for (int i = 0; i < diffDelta.size(); ++i)
//	{
//		printf(" %.1lf", predictedDiffDelta(i));
//	}
//	printf("\n");
//	printf("Act ");
//	for (int i = 0; i < diffDelta.size(); ++i)
//	{
//		printf(" %.1lf", diffDelta(i));
//	}
//	printf("\n");
		system.zVecInitial = tmp + h * 1/2 * k2;
		VectorXd k3 = getAdjustVector(system, matShuffle, kernelMatrix);
		if (k3.size() == 0) { printf("Reset!\n"); system = System(system.num); continue; }
		system.zVecInitial = tmp + h * 1 * k3;
		VectorXd k4 = getAdjustVector(system, matShuffle, kernelMatrix);
		if (k4.size() == 0) { printf("Reset!\n"); system = System(system.num); continue; }

		system.zVecInitial = tmp + h * (1/6. * k1 + 1/3. * k2 + 1/3. * k3 + 1/6. * k4);
	}

	/*System temp = system;
	evolve(temp, false);
	for (int i = 0; i < temp.params; ++i)
	{
		for (int j = 0; j < temp.params; ++j)
		{
			printf(" %.1lf", temp.aMatTotal(i, j));
		}
		printf("\n");
	}
	printf("\n");
	VectorXd diff = temp.zVec - temp.zVecInitial;
	VectorXd adjust;
	temp.aMatTotal.lu().solve(diff, &adjust);
	for (int i = 0; i < temp.params; ++i)
	{
		printf(" %.1lf", temp.zVecInitial(i));
	}
	printf("\n");
	for (int i = 0; i < temp.params; ++i)
	{
		printf(" %.1lf", temp.zVec(i));
	}
	printf("\n");
	for (int i = 0; i < temp.params; ++i)
	{
		printf(" %.1lf", diff(i));
	}
	printf("\n");
	for (int i = 0; i < temp.params; ++i)
	{
		printf(" %.1lf", adjust(i));
	}
	printf("\n");*/

	/*system.zVec = system.zVecInitial;
	for (int d = 0; d < D; ++d)
	{
		double average = 0;
		for (int i = 0; i < system.num; ++i)
		{
			average += system.Q(i, d)/system.num;
		}
		for (int i = 0; i < system.num; ++i)
		{
			system.Q(i, d) -= average;
		}
	}
	system.zVecInitial = system.zVec;*/

	Image beforeAfter(imageSize, imageSize, 0 - border, 0 - border, 1 + border, 1 + border);
	system.render(beforeAfter);

	evolve(system, true);

	system.render(beforeAfter);
	beforeAfter.save("BeforeAfter");

	/*for (int i = 0; i < system.params; ++i)
	{
		for (int j = 0; j < system.params; ++j)
		{
			printf(" %.12lf", system.aMatTotal(i, j));
		}
		printf("\n");
	}
	printf("\n");*/
	for (int i = 0; i < system.zVecInitial.size(); ++i)
	{
		printf("  %.2lf", system.zVecInitial(i));
	}
	printf("\n");
	for (int i = 0; i < system.zVec.size(); ++i)
	{
		printf("  %.2lf", system.zVec(i));
	}
	printf("\n");
}
