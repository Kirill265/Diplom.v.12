#include <iostream> // для организации ввода-вывода
#include <string>   // для организации работы со строками 
#include <fstream>  // для считывания и записи данных из файла и в файл
#include <iomanip>  // для setprecision()
#include <cmath>    // математические выражения

using namespace std;//подключение пространства стандартных имен

double F[11] = { 0.0 }, Y0[30], Y1[30], X0[30], data1[14500][40], data2[14500][40], A[3][3], B[3][3], Matrix[3][3], t = 0, shag;
double aa[5], b[5], c[2], t_zad, X_t, t_border, tangaj_1st, ryskanie_1st, kren_1st, delta_teta, delta_psi;
double promah, epsilon, Promah[14500], Epsilon[14500], y_tg, x_tg, VxTG1, VyTG1, VxTG2, VyTG2, alfaTG, AlphaTG[14500], mu_new;

const double pi = 3.141592;
unsigned int choise = 0;  // флаг для меню
unsigned int rocket = 0;  // флаг для выбора ракеты
unsigned int target = 0;  // флаг для выбора цели
unsigned int choise2 = 0; // флаг для меню
unsigned int exper = 0;	  // флаг для расчетов с при прогнозе
unsigned int prognoz = 1;	  // флаг для прогноза
double F0[11] = { 0.0 };    // { фи, хи, альфа, бетта, del_v, del_n, del_el, радиус визирования, перегрузка, угол наклона траектории, угол пути}

double V0 = 350;					// скорость на входе из шахты
double tangaj0 = 30.0 * pi / 180.0; // тангаж
double ryskanie0 = 0.0;				// рыскание
double kren0 = 0.0;					// крен
double mass0 = 980.0;				// масса
double alf_Tg[3] = { 0.0 }, bet_Tg[3] = { 0.0 }; // угол альфа и бетта

								  //перевод из углов элера в кватернионы
double Ro = cos(ryskanie0 / 2)*cos(tangaj0 / 2)*cos(kren0 / 2) - sin(ryskanie0 / 2)*sin(tangaj0 / 2)*sin(kren0 / 2);
double Lambda = sin(ryskanie0 / 2)*sin(tangaj0 / 2)*cos(kren0 / 2) + cos(ryskanie0 / 2)*cos(tangaj0 / 2)*sin(kren0 / 2);
double Mu = sin(ryskanie0 / 2)*cos(tangaj0 / 2)*cos(kren0 / 2) + cos(ryskanie0 / 2)*sin(tangaj0 / 2)*sin(kren0 / 2);
double Nu = cos(ryskanie0 / 2)*sin(tangaj0 / 2)*cos(kren0 / 2) - sin(ryskanie0 / 2)*cos(tangaj0 / 2)*sin(kren0 / 2);

// для целей и ракет
//                  0             1                     2            3    4     5   6    7          8      9      10  11      12  13   14     15		 16			 17       18       19         20         21        22         23        24        25         26         27        28          29
//               {  ?,           Vxg,                  Vyg,         Vzg,  Xg,  Yg, Zg, Omega x, Omega y, Omega z, Ro, Lambda, Mu, Nu, kren,	 ryskanie,	tangaj,		масса, x цели 1, y цели 1, z цели 1, Vx цели 1, Vy цели 1, Vz цели 1, x цели 2, y цели 2, z цели 2 , Vx цели 2, Vy цели 2, Vz цели 2 }
double X00[30] = { 0.0, V0 * cos(tangaj0),	V0 * sin(tangaj0),		0.0, 0.0, 0.0, 0.0, 0.0,      0.0,     0.0,   Ro, Lambda, Mu, Nu, kren0, ryskanie0, tangaj0,	mass0,  23523.0,	5000.0,	 500.0,    -200.00,		0.0,    0.0,     32932.0,	8000.0,    0.0,     -250.0,       0.0,       0.0 };
 
//для прогноза полета к 1-ой и 2-ой целям:
double Xf[30] = { 0.0 }, Xs[30] = { 0.0 };	//параметры ракеты и траектории
double tf = 0, ts = 0;						//время полета
double F_M[11] = { 0.0 };					//сохранение в памяти параметра F
double dr_dt1[14500] = { 0.0 }, dr_dt2[14500] = { 0.0 };	//производные по времени радиусов визирования до цели 1 и 2 сосответсвенно
double var_H, var_L;						//изменение паарметров цели наведения
double tf_m, ts_m;							//момент, когда перегрузка больше нормы
double n_Ya_toF, n_Ya_toS;					//максимальная перегрузка при полете к первой и второй целям среди всех точек траектории
double n_Ya_1[14500] = { 0.0 }, n_Ya_2[14500] = { 0.0 };	//максимальная перегрузка при полете к первой и второй целям в каждой точке траектории
double r_viz_f[14500] = { 0.0 }, r_viz_s[14500] = { 0.0 };	//конечный радиус визирования при полете к первой и второй целям
double r_viz_l1 = 0, r_viz_l2 = 0;			//последний раситанный конечный радиус визирования при полете к первой и второй целям


double dRo = 0.0;
double dLambda = 0.0;
double dMu = 0.0;
double dNu = 0.0;

double Tangaj = 0.0;
double Kren = 0.0;
double Ryskanie = 0.0;

//интегрирование параметров Родриго-Гамильтона
void IntRG(const double Ro, const double Lambda, const double Mu, const double Nu, const double OmegaX, const double OmegaY, const double OmegaZ)
{
	dRo = -(OmegaX*Lambda + OmegaY*Mu + OmegaZ*Nu) / 2;
	dLambda = (OmegaX*Ro - OmegaY*Nu + OmegaZ*Mu) / 2;
	dMu = (OmegaX*Nu + OmegaY*Ro - OmegaZ*Lambda) / 2;
	dNu = (-OmegaX*Mu + OmegaY*Lambda + OmegaZ*Ro) / 2;
}

//первод из кватернионов в углы Эйлера
void Ugl(const double Ro, const double Lambda, const double Mu, const double Nu)
{
	Ryskanie = atan(2 * (Ro*Mu - Lambda*Nu) / (pow(Ro, 2) + pow(Lambda, 2) - pow(Mu, 2) - pow(Nu, 2)));
	Kren = atan(2 * (Ro*Lambda - Mu*Nu) / (pow(Ro, 2) - pow(Lambda, 2) + pow(Mu, 2) - pow(Nu, 2)));
//	Ryskanie = atan2(2 * (Ro*Mu - Lambda*Nu), pow(Ro, 2) + pow(Lambda, 2) - pow(Mu, 2) - pow(Nu, 2));
//	Kren = atan2(2 * (Ro*Lambda - Mu*Nu), pow(Ro, 2) - pow(Lambda, 2) + pow(Mu, 2) - pow(Nu, 2));
	Tangaj = asin(2 * (Ro*Nu + Mu*Lambda));
}

//поиск минимального из чисел
double min(double x1, double x2)
{
	if (x1 <= x2)
		return x1;
	else //x2 < x1
		return x2;
}
// поиск минимального из чисел
double max(double x1, double x2)
{
	if (x1 >= x2)
		return x1;
	else //x2 > x1
		return x2;
}
double closer(double x1, double x2, double p1, double p2)
{
	if (x1 <= x2)
		return p1;
	else //x2 < x1
		return p2;
}
double further(double x1, double x2, double p1, double p2)
{
	if (x1 >= x2)
		return p1;
	else //x2 < x1
		return p2;
}


// прототипы функций
double sk_zv(double height);        //Функция расчета скорсоти звука в атмосфере
double pl(double height);			//Функция расчета плотности атмосферы
double cx(double M, double alfa_);  //Функции нахождения АДХ атмосферы (включая производные по углам альфа, бетта, дельта)
double cy_a(double M, double alfa_);
double cz_b(double M, double alfa_);
double cy_d(double M, double alfa_);
double cz_d(double M, double alfa_);
double mx_wx(double M, double alfa_);
double mx_d(double M, double alfa_);
double mz_wz(double M, double alfa_);
double mz_a(double M, double alfa_);
double mz_d(double M, double alfa_);
double my_wy(double M, double alfa_);
double my_b(double M, double alfa_);
double my_d(double M, double alfa_);
double theta_sr(double theta_closer, double V_closer, double xn_closer, double theta_further, double V_further, double yn_further);

double viz_fi(double Y, double r_viz, double Yc);        //Функция определния наклона линии визирования в вертикальной плоскости
double viz_hi(double X, double Xc, double Z, double Zc); //Функция определния наклона линии визирования в горизонтальной плоскости

double r_vizV(double X, double Y, double Z, double Xc, double Yc, double Zc); //Функция определния радиуса линии визирования

void Eyler(double tau, double X[30]);
void RG_sv_st(double Ro, double Lambda, double Mu, double Nu);
void A_sv_st(double t, double r, double k);
void B_sf_st(double fi, double hi);
void M_sv_sk(double alfa, double betta);
void Miss(double Y_LA, double X_LA, double theta_la, double V_la, double Y_TG, double X_TG, double theta_tg, double V_tg);
void Prognoz(double y_new_tg, double x_new_tg, double Vy_new_tg, double Vx_new_tg, double Y_LA, double X_LA, double theta_la);

int main()
{
	setlocale(LC_ALL, "Russian");

	//присвоение начальных данных
	int sh = 0;
	int i_f = 0, i_s = 0;
	// Вывод параметров в виде таблицы
	string X_names[38] = {
		"Vxg                ", "Vyg                ", "Vzg                ",
		"Xg                 ", "Yg                 ", "Zg                 ",
		"Omega x            ", "Omega y            ", "Omega z            ",
		"Ro                 ", "Lambda             ", "Mu                 ",
		"Nu                 ", "крен               ", "рыскание           ",
		"тангаж             ", "масса              ", "x цели 1           ",
		"y цели 1           ", "z цели 1           ", "Vx цели 1          ",
		"Vy цели 1          ", "Vz цели 1          ", "x цели 2           ",
		"y цели 2           ", "z цели 2           ", "Vx цели 2          ",
		"Vy цели 2          ", "Vz цели 2          ", "фи                 ",
		"хи                 ", "альфа              ", "бетта              ",
		"del_v              ", "del_n              ", "del_el             ",
		"радиус визирования ", "перегрузка         "
	};
	shag = 0.01;

	do
	{
		t = 0;
		sh = 0;
		var_H = 0;
		var_L = 0;
		for (int i = 0; i < 30; i++)
		{
			X0[i] = X00[i];
			if (i >= 0 && i < 11)
			{
				F[i] = F0[i];
			}
		}
		do
		{   // выбор ракеты
			cout << "Выберите одну из двух ракет. Нажмите 1 или 2.\n";
			cin >> rocket;
			if (rocket != 1 && rocket != 2)
				cout << "Был произведен некорректный ввод данных!\n";
		} while (rocket != 1 && rocket != 2);

		if (rocket == 1)
		{
			tangaj_1st = tangaj0;
			ryskanie_1st = ryskanie0;
		}
		else // rocket == 2
		{
			do
			{   // выбор цели
				cout << "Выберите одну из двух целей. Нажмите 1 или 2.\n";
				cin >> target;
				if (target != 1 && target != 2 && target != 0)
					cout << "Был произведен некорректный ввод данных!\n";
			} while (target != 1 && target != 2 && target != 0);

//			cout << "Введите интервал запуска: "; //ввод интервала между пусками 1-2 ракет
//			cin >> t_zad;
			t_zad = 2;

			//отличающиеся начальные параметры ракет
			tangaj0 = 30 * pi / 180;
			X0[16] = tangaj0;						//начальный угол тангажа для второй ракеты
			//перерасчёт для новых углов
			X0[1] = V0 * cos(tangaj0);
			X0[2] = V0 * sin(tangaj0);
			Ro = cos(ryskanie0 / 2)*cos(tangaj0 / 2)*cos(kren0 / 2) - sin(ryskanie0 / 2)*sin(tangaj0 / 2)*sin(kren0 / 2);
			Lambda = sin(ryskanie0 / 2)*sin(tangaj0 / 2)*cos(kren0 / 2) + cos(ryskanie0 / 2)*cos(tangaj0 / 2)*sin(kren0 / 2);
			Mu = sin(ryskanie0 / 2)*cos(tangaj0 / 2)*cos(kren0 / 2) + cos(ryskanie0 / 2)*sin(tangaj0 / 2)*sin(kren0 / 2);
			Nu = cos(ryskanie0 / 2)*sin(tangaj0 / 2)*cos(kren0 / 2) - sin(ryskanie0 / 2)*cos(tangaj0 / 2)*sin(kren0 / 2);
			X0[10] = Ro;
			X0[11] = Lambda;
			X0[12] = Mu;
			X0[13] = Nu;
			tangaj_1st = X0[16];
			ryskanie_1st = X0[15];
			kren_1st = X0[14];
			X0[18] += X0[21] * t_zad;
			X0[19] += X0[22] * t_zad;
			X0[20] += X0[23] * t_zad;
			X0[24] += X0[27] * t_zad;
			X0[25] += X0[28] * t_zad;
			X0[26] += X0[29] * t_zad;
		}

		do
		{
/*			if (F[7] > 50 && rocket == 2 && prognoz == 1)
			{
				exper = 1;
				if (n_Ya_toF <= 15 && r_viz_l1 <= 50)
			//	if (true)
				{
					tf = t;
					for (int i = 0; i < 30; i++)
					{
						Xf[i] = X0[i];
						if (i >= 0 && i < 11)
						{
							F_M[i] = F[i];
						}
					}
					target = 1;
					i_f = -1;
					do
					{
						i_f += 1;
						Eyler(tf, Xf);
						for (int i = 1; i < 30; i++)
						{
							Xf[i] += shag*Y0[i];
						}
						tf += shag;
						if (F[8] > n_Ya_1[sh])
							n_Ya_1[sh] = F[8];
						if (i_f == 1)
							dr_dt1[sh] -= F[7];
						r_viz_l1 = F[7];
					} while ((F[7] >= 50) && (Xf[24] >= 1000) && (Xf[5] >= 0) && (Xf[4] <= Xf[18]) && (Xf[16] <= pi) && (tf <= 48));
					dr_dt1[sh] += r_viz_l1;
					dr_dt1[sh] /= (i_f - 1);
					r_viz_f[sh] = r_viz_l1;
					if (n_Ya_1[sh] > n_Ya_toF)
						n_Ya_toF = n_Ya_1[sh];

					tf_m = t;

					for (int i = 0; i < 11; i++)
					{
					F[i] = F_M[i];
					}
				}

				if (n_Ya_toS <= 15 && r_viz_l2 <= 50)
			//	if (true)
				{
					ts = t;
					for (int i = 0; i < 30; i++)
					{
						Xs[i] = X0[i];
						if (i >= 0 && i < 11)
						{
							F_M[i] = F[i];
						}
					}
					target = 2;
					i_s = -1;
					do
					{
						i_s += 1;
						Eyler(ts, Xs);
						for (int i = 1; i < 30; i++)
						{
							Xs[i] += shag*Y0[i];
						}
						ts += shag;
						if (F[8] > n_Ya_2[sh])
							n_Ya_2[sh] = F[8];
						if (i_s == 1)
							dr_dt2[sh] -= F[7];
						r_viz_l2 = F[7];
					} while ((F[7] >= 50) && (Xs[24] >= 1000) && (Xs[5] >= 0) && (Xs[4] <= Xs[24]) && (Xs[16] <= pi) && (ts <= 48));
					dr_dt2[sh] += r_viz_l2;
					dr_dt2[sh] /= (i_s - 1);
					r_viz_s[sh] = r_viz_l2;
					if (n_Ya_2[sh] > n_Ya_toS)
						n_Ya_toS = n_Ya_2[sh];

					ts_m = t;

					for (int i = 0; i < 11; i++)
					{
					F[i] = F_M[i];
					}
				}
				if (dr_dt1[sh] == 0)
				{
					if (dr_dt2[sh] == 0) //все в норме, тот же курс
					{
						cout << "\n\tAHTUNG!!!\nперегрузка превышает норму по двум направлениям:\nк 1 цели:  " << tf_m << "\tко 2 цели:  " << ts_m << "\n";
					}
					else //ко второй ракете
					{
						//cout << "\nцель: вторая ракета:\nмаксимальная перегрузка к 1 цели:  " << n_Ya_toF << "\n";
						target = 2;
						prognoz = 0;
					}
				}
				else
				{
					if (dr_dt2[sh] == 0) //к первой ракете
					{
						//cout << "\nцель: первая ракета:\nмаксимальная перегрузка ко 2 цели:  " << n_Ya_toS << "\n";
						target = 1;
						prognoz = 0;
					}
					else //корректируем параметры наведения
					{
						if (dr_dt2[sh] > dr_dt1[sh]) //увеливаем угол наклона радиуса визирования (смещаем точку наведения ближе к нам или выше)
						{
						//	var_H += ;
						//	var_L += ;
							target = 0;
						}
						else if (dr_dt2[sh] < dr_dt1[sh]) //уменьшаем угол наклона радиуса визирования (смещаем точку наведения дальше от нас или ниже)
						{
						//	var_H -= ;
						//	var_L -= ;
							target = 0;
						}
						else //dr_dt2[sh] == dr_dt1[sh]   //все отлично!
						{
							target = 0;
						}
					}
				}
				exper = 0;
			}
*/			
			//интегрирование
			Eyler(t, X0); //получение правых частей
			for (int i = 1; i < 30; i++)
			{
				X0[i] += shag*Y0[i];
			}
			X0[21] = VxTG1;
			X0[22] = VyTG1;
			X0[27] = VxTG2;
			X0[28] = VyTG2;
			t = t + shag;
			//запись текущих параметров в матрицу
			for (int i = 0; i < 30; i++)
			{
				data2[sh][i] = X0[i];
				if (i >= 0 && i<11)
					data2[sh][i + 30] = F[i];
			}
			AlphaTG[sh] = alfaTG;
			Promah[sh] = promah;
			Epsilon[sh] = epsilon;
			if (rocket == 1)
			{
				X_t = X0[18]; //положение цели
//				t_border = t; //время подрыва первой ракеты
			}
			else if (target == 1) // rocket == 2
			{
				X_t = X0[18];
			}
			else // rocket == 2 && target == 2
			{
				X_t = X0[24];
			}
			if (F[7] > 10000)
			{
				t_border = t; //время подрыва первой ракеты
			}
			sh = sh + 1;
		} while (((F[7] >= 50 && rocket == 2) || (F[7] >= 50 && rocket == 1)) && X_t >= 1000 && X0[5] >= 0 && X0[4] <= X_t && X0[16] <= pi && t <= 148);


		cout << "Конечные параметры ракеты и целей : \n";
		//вывод конечных параметров ракеты и целей
		for (int i = 1; i < 30; i++)
		{
			cout << X_names[i - 1] << " " << X0[i] << "\n";
		}
		cout << "r_viz               " << F[7] << "\n";
		cout << "t                   " << t << "\n";
		cout << "alf_Tg              " << alf_Tg[0] << "\n";
		cout << "bet_Tg              " << bet_Tg[0]  << "\n";
		cout << "bet_Tg              " << bet_Tg[0] << "\n";
//		cout << "t_border            " << t_border << "\n";
		cout << "n_Ya1               " << n_Ya_toF << "\t" << "r_V1     " << r_viz_l1 << "\t" << "t    " << tf_m << "\n";
		cout << "n_Ya2               " << n_Ya_toS << "\t" << "r_V2     " << r_viz_l2 << "\t" << "t    " << ts_m << "\n";
		do
		{
			cout << "Нажмите 1, если хотите повторить расчет, иначе - 0.\n";
			cin >> choise;
			if (choise != 0 && choise != 1)
				cout << "Был произведен некорректный ввод данных!\n";
		} while (choise != 1 && choise != 0);

		if (choise == 0)
		{
			do
			{
				cout << "Хотите ли записать результаты расчетов в файл, если да, то нажмите 1, иначе - 0.\n";
				cin >> choise2;
				if (choise2 != 0 && choise2 != 1)
					cout << "Был произведен некорректный ввод данных!\n";
			} while (choise2 != 1 && choise2 != 0);

			//запись параметров в файл
			ofstream ffout("3D_02.txt");
			if (choise2 == 1)
			{
//				ffout << "t, c \t" << "V, м/с \t" << "tetta, град\t" << "X, м\t" << "Y, м\t" << "Z, м\t" << "omegaZ, рад/c\t"
//					<< "m, kg\t" << "Xc_1, м\t" << "Yc_1, м\t" << "Zc_1, м\t" << "Xc_2, м\t" << "Yc_2, м\t" << "Zc_2, м\t" << "n_Ya, -\t"
//					<< "fi, град\t" << "hi, град\t" << "r_viz м\t" << "alfa, град\t" << "betta, град\t" << "delta_v, град\t" << "delta_n, град\n";
				for (int i = 0; i < sh; i += 10)
				{
					Ugl(data2[i][10], data2[i][11], data2[i][12], data2[i][13]);

					ffout /*<< fixed << setprecision(2) */ << i*shag << "\t";
					ffout /*<< fixed << setprecision(4)*/ << sqrt(pow(data2[i][1], 2) + pow(data2[i][2], 2) + pow(data2[i][3], 2)) << "\t";
					ffout << fixed << setprecision(4) << Tangaj * 180 / 3.141592 << "\t"; 
					ffout /*<< fixed << setprecision(4)*/ << data2[i][4] << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][5] << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][6] << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][9] << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][17] << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][18] << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][19] << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][20] << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][24] << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][25] << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][26] << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][38] << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][30] * 180 / 3.141592 << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][31] * 180 / 3.141592 << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][37] << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][32] * 180 / 3.141592 << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][33] * 180 / 3.141592 << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][34] * 180 / 3.141592 << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][35] * 180 / 3.141592 << "\t";
					ffout << dr_dt1[i] << "\t";
					ffout << dr_dt2[i] << "\t";
					ffout << r_viz_f[i] << "\t";
					ffout << r_viz_s[i] << "\t";
					ffout << n_Ya_1[i] << "\t";
					ffout << n_Ya_2[i] << "\t";
					ffout << Promah[i] << "\t";
					ffout << Epsilon[i] * 180 / 3.141592 << "\t";
					ffout << AlphaTG[i] << "\t";
					ffout << fixed << setprecision(4) << Ryskanie * 180 / 3.141592 << endl;


				}
				sh = sh - 1;
				Ugl(data2[sh][10], data2[sh][11], data2[sh][12], data2[sh][13]);
				ffout /*<< fixed << setprecision(2) */ << sh*shag << "\t";
				ffout /*<< fixed << setprecision(4)*/ << sqrt(pow(data2[sh][1], 2) + pow(data2[sh][2], 2) + pow(data2[sh][3], 2)) << "\t";
				ffout << fixed << setprecision(4) << Tangaj * 180 / 3.141592 << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][4] << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][5] << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][6] << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][9] << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][17] << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][18] << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][19] << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][20] << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][24] << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][25] << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][26] << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][38] << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][30] * 180 / 3.141592 << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][31] * 180 / 3.141592 << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][37] << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][32] * 180 / 3.141592 << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][33] * 180 / 3.141592 << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][34] * 180 / 3.141592 << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][35] * 180 / 3.141592 << "\t";
				ffout << dr_dt1[sh] << "\t";
				ffout << dr_dt2[sh] << "\t";
				ffout << r_viz_f[sh] << "\t";
				ffout << r_viz_s[sh] << "\t";
				ffout << n_Ya_1[sh] << "\t";
				ffout << n_Ya_2[sh] << "\t";
				ffout << Promah[sh] << "\t";
				ffout << Epsilon[sh] * 180 / 3.141592 << "\t";
				ffout << AlphaTG[sh] << "\t";
				ffout << fixed << setprecision(4) << Ryskanie * 180 / 3.141592 << endl;
			}
			do
			{
				cout << "\nПроизвести еще какие-то расчеты?\nДля выхода из программы нажмите 0, для проведения расчетов нажмите 1.\n";
				cin >> choise;
				if (choise != 0 && choise != 1)
					cout << "Был произведен некорректный ввод данных!\n";
			} while (choise != 0 && choise != 1);
		}
	} while (choise == 1);
	return 0;
}
//функция нахождения правых частей ДУ
void Eyler(double tau, double X[30])
{
	double d = 0.5;
	double Ix = 240;
	double Iy = 1500;//2239;//2239; //640
	double Iz = 1500;//2239; //640
	double m = X[17];
	double g = 9.80665;
	double a = sk_zv(X[5]);
	double ro = pl(X[5]);
	double l = 6.3;
	double t_1st = 2.00;					//продолжительность программного участка траектории
//	double t_gambit = 5;				//момент манёвра
	Y0[17] = -18.174 + 0.375*tau;		//закон изменения массового расхода
	if (tau >= 48)						//условие конца горения топлива
		Y0[17] = 0;

	double P = 3230.0 * (-Y0[17]);

	double deltaAlf_Tg[3] = { 0, 0, 0 };
	double deltaBet_Tg[3] = { 0, 0, 0 };
	double Vx_newTg[3] = { 0, 0, 0 }, Vy_newTg[3] = { 0, 0, 0 }, Vz_newTg[3] = { 0, 0, 0 };
	double V_Tg[3] = { 0, 0, 0 }, Theta_Tg[3] = { 0, 0, 0 };

	double r_vizir1 = r_vizV(X[4], X[5], X[6], X[18], X[19], X[20]);

		//матрица перевода из скоростной СК в связанную
		M_sv_sk(deltaAlf_Tg[0], deltaBet_Tg[0]);
		//изменение наклона траектории цели
		Vx_newTg[0] = Matrix[0][0] * X[21] + Matrix[0][1] * X[22] + Matrix[0][2] * X[23];
		Vy_newTg[0] = Matrix[1][0] * X[21] + Matrix[1][1] * X[22] + Matrix[1][2] * X[23];
		Vz_newTg[0] = Matrix[2][0] * X[21] + Matrix[2][1] * X[22] + Matrix[2][2] * X[23];

		//матрица перевода из скоростной СК в связанную
		M_sv_sk(deltaAlf_Tg[1], deltaBet_Tg[1]);
		//изменение наклона траектории цели
		Vx_newTg[1] = Matrix[0][0] * X[27] + Matrix[0][1] * X[28] + Matrix[0][2] * X[29];
		Vy_newTg[1] = Matrix[1][0] * X[27] + Matrix[1][1] * X[28] + Matrix[1][2] * X[29];
		Vz_newTg[1] = Matrix[2][0] * X[27] + Matrix[2][1] * X[28] + Matrix[2][2] * X[29];

/*
		alf_Tg[1] += deltaAlf_Tg[1];
		bet_Tg[1] += deltaBet_Tg[1];

	}
	if ((rocket == 2) && (tau < t_border - t_zad)) //условие: для 2-ой ракеты (когда первая ракета еще не взорвалась)
	{
		//берем данные о целях из памяти 
		int t_otn = int((tau + t_zad) / shag);
		for (int i = 18; i < 21; i++)
			X[i] = data2[t_otn - 1][i];
		for (int i = 21; i < 24; i++)
			X[i] = data2[t_otn][i];
		for (int i = 24; i < 27; i++)
			X[i] = data2[t_otn - 1][i];
		for (int i = 27; i < 30; i++)
			X[i] = data2[t_otn][i];
	}
*/

	Ugl(X[10], X[11], X[12], X[13]);

	//нормирование кватернионов
	double L1 = sqrt(pow(X[10], 2) + pow(X[11], 2) + pow(X[12], 2) + pow(X[13], 2));
	for (int i = 10; i < 14; i++)
		X[i] /= L1;

	double tangaj = Tangaj;
	double ryskanie = Ryskanie;
	double kren = Kren;

	double Vxg = X[1];
	double Vyg = X[2];
	double Vzg = X[3];

	double V = sqrt(pow(Vxg, 2) + pow(Vyg, 2) + pow(Vzg, 2));
	double q = ro*pow(V, 2) / 2;
	double S = pi*pow(d, 2) / 4;
	double G = -m *g;
	double M = V / a;

	//матрица переовода из стартовой СК в связанную
	RG_sv_st(X[10], X[11], X[12], X[13]);

	//перевод скоростей ракеты из стартовой СК в связанную
	double Vx = A[0][0] * Vxg + A[0][1] * Vyg + A[0][2] * Vzg;
	double Vy = A[1][0] * Vxg + A[1][1] * Vyg + A[1][2] * Vzg;
	double Vz = A[2][0] * Vxg + A[2][1] * Vyg + A[2][2] * Vzg;

	double alfa = -atan2(Vy, Vx);
	double betta = asin(Vz / V);
	double Q = asin(Vyg / V);						//угол наклона траеткории
	double Put = atan2(-Vzg, Vxg);					//угол пути
	double alfa1 = sqrt(pow(alfa, 2) + pow(betta, 2));

	//шаг изменения углов наклона траектории
	delta_teta = Q - F[9];
	delta_psi = Put - F[10];
//	cout << delta_teta << "\n";
	F[2] = alfa;
	F[3] = betta;
	F[9] = Q;
	F[10] = Put;

	//расчет АДХ
	double Cx = cx(M, alfa1);
	double Cy_a = cy_a(M, alfa); // * 1.5
	double Cz_b = cz_b(M, betta);// / 5;

	double Cydelta = cy_d(M, alfa) * 1.5;
	double Czdelta = cz_d(M, betta);// / 5;

	double mx_w = mx_wx(M, alfa1)*q*S*l;
	double my_w = my_wy(M, betta)*q*S*l; // *5
	double mz_w = mz_wz(M, alfa)*q*S*l;

	double my_betta = my_b(M, betta)*q*S*l;
	double mz_alfa = mz_a(M, alfa)*q*S*l;

	double Mxdel = mx_d(M, alfa1)*q*S*l;
	double Mydel = my_d(M, betta)*q*S*l * 2;
	double Mzdel = mz_d(M, alfa)*q*S*l * 0.5;

	//получение динамических парамтеров наведения и стабилизации в момент полета
	
	aa[1] = -mz_w*l / (Iz*V);
	aa[2] = -mz_alfa / Iz;
	aa[3] = -Mzdel / Iz;
	aa[4] = (P + Cy_a*q*S) / (m*V);
	aa[5] = Cydelta*q*S / (m*V);

	b[1] = -my_w*l / (Iy*V);
	b[2] = -my_betta / Iy;
	b[3] = -Mydel / Iy;
	b[4] = (-P + Cz_b*q*S) / (m*V);
	b[5] = -Czdelta*q*S / (m*V);
//	b[4] = (P + Cz_b*q*S) / (m*V);
//	b[5] = Czdelta*q*S / (m*V);

	c[1] = -mx_w*l / (Ix*V);
	c[3] = -Mxdel / Ix;

	double ksiCC_v = 0.35;
	double KCC_v = 0.95;
	double ksiCC_n = 0.35;
	double KCC_n = 0.95;
	double ksiCC_el = 0.35;
	double TCC_el = 0.01;

	double K_v = (aa[2] * aa[5] - aa[3] * aa[4]) / (aa[2] + aa[1] * aa[4]);
	double T1_v = -aa[3] / (aa[3] * aa[4] - aa[2] * aa[5]);
	double T_v = 1 / sqrt(aa[2] + aa[1] * aa[4]);
	double ksi_v = (aa[1] + aa[4]) / (2 * sqrt(aa[2] + aa[1] * aa[4]));

	double K_n = (b[2] * b[5] - b[3] * b[4]) / (b[2] + b[1] * b[4]);
	double T1_n = -b[3] / (b[3] * b[4] - b[2] * b[5]);
	double T_n = 1 / sqrt(b[2] + b[1] * b[4]);
	double ksi_n = (b[1] + b[4]) / (2 * sqrt(b[2] + b[1] * b[4]));

	double K_el = c[3] / c[1];
	double T_el = 1 / c[1];

	double K2_v = -2 * T_v*(ksi_v*T1_v - pow(ksiCC_v, 2)*T_v - sqrt(pow(T_v, 2)*pow(ksiCC_v, 4) - 2 * ksi_v*pow(ksiCC_v, 2)*T_v*T1_v + pow(T1_v*ksiCC_v, 2))) / (K_v*pow(T1_v, 2));
	double K1_v = KCC_v*(1 + K_v*K2_v) / (K_v*pow(T1_v, 2));

	double K2_n = -2 * T_n*(ksi_n*T1_n - pow(ksiCC_n, 2)*T_n - sqrt(pow(T_n, 2)*pow(ksiCC_n, 4) - 2 * ksi_n*pow(ksiCC_n, 2)*T_n*T1_n + pow(T1_n*ksiCC_n, 2))) / (K_n*pow(T1_n, 2));
	double K1_n = KCC_n*(1 + K_n*K2_n) / (K_n*pow(T1_n, 2));

	double K1_el = (2 * ksiCC_el*T_el - TCC_el) / (K_el*TCC_el);
	double K2_el = T_el / (K_el*pow(TCC_el, 2));

	double K_fi = 50.55;
	double Kk_v = 30;

	double K_hi = +0.5;
	double Kk_n = 1;

	double Xc, Yc, Zc;
	double Vxc, Vyc, Vzc;

	//параметры цели для первой ракеты
	if (rocket == 1)
	{
		if (tau < t_1st) //на участке выведения
		{
			/*
			Xc = X[18];
			Yc = X[19];
			Zc = X[20];
			Vxc = X[21];
			Vyc = X[22];
			Vzc = X[23];
			*/
			Xc = 30000*cos(tangaj_1st)*cos(ryskanie_1st);
			Yc = 30000*sin(tangaj_1st);
			Zc = 30000*cos(tangaj_1st)*sin(ryskanie_1st);
			Vxc = 0;
			Vyc = 0;
			Vzc = 0;
			
		}
		else if (r_vizir1 >= 6000) //на втором участке траектории
		{
			Xc = X[18];
			Yc = X[19];
			Zc = X[20];
			Vxc = X[21];
			Vyc = X[22];
			Vzc = X[23];
		}
		else //на участке минимакс
		{
			Xc = X[18];
			Yc = X[19];
			Zc = X[20];
			Vxc = X[21];
			Vyc = X[22];
			Vzc = X[23];
		}
	}
	// задание параметров средней точки наведения до выбора конкретной цели
	double Xc_12 = min(X[18], X[24]); // -var_L;
//	double Yc_12 = max(X[19], X[25]) + var_H;

//	double Xc_12 = 0.5*(X[18] + X[24]);
	double Yc_12 = 0.5*(X[19] + X[25]);
	double Zc_12 = 0.5*(X[20] + X[26]);
	
	double Vxc_12 = closer(X[18], X[24], X[21], X[27]);
	double Vyc_12 = 0.5*X[22] + 0.5*X[28];
	double Vzc_12 = 0.5*X[23] + 0.5*X[29];

	double R_viz_12 = r_vizV(X[4], X[5], X[6], Xc_12, Yc_12, Zc_12);
	//выбор цели
	if (rocket == 2)
	{
		if (tau < t_1st && exper == 0) //на программном участке траектории
		{
			Xc = 30000*cos(tangaj_1st)*cos(ryskanie_1st);
			Yc = 30000*sin(tangaj_1st);
			Zc = 30000*cos(tangaj_1st)*sin(ryskanie_1st);
			Vxc = 0;
			Vyc = 0;
			Vzc = 0;
		}
//		else if (target == 0)
		else if (R_viz_12 >= 10000)
		//		else if (t <= t_border - t_zad) //на втором участке траектории
		{
			/*
			Xc = X[18];
			Yc = X[19];
			Zc = X[20];
			Vxc = X[21];
			Vyc = X[22];
			Vzc = X[23];
			*/
			Xc = Xc_12;
			Yc = Yc_12;
			Zc = Zc_12;
			Vxc = Vxc_12;
			Vyc = Vyc_12;
			Vzc = Vzc_12;
			
		}
		else //на участке траектории с конкретной целью
			 //			if (R_viz_12 >= 30) 
		{
			if (target == 1)
			{
				Xc = 15000; // X[18];
				Yc = 6000;// X[19];
				Zc = X[20];
				Vxc = X[21];
				Vyc = X[22];
				Vzc = X[23];
			}
			if (target == 2)
			{
				Xc = 19000; // X[24];
				Yc = 8000; //X[25];
				Zc = X[26];
				Vxc =  X[27];
				Vyc =   X[28];
				Vzc =   X[29];
			}
		}
	}
	F[7] = r_vizV(X[4], X[5], X[6], Xc, Yc, Zc);
	
	//////////////////////////////////
		double y_new_tg1 = X[19];
		double x_new_tg1 = X[18];
		double y_new_tg2 = X[25];
		double x_new_tg2 = X[24];
		double deltaAlf_pr1 = 0;
		if (r_vizV(X[4], X[5], X[6], X[18], X[19], X[20]) < 15500)
		{
			Miss(X[5], X[4], Q, V, X[19], X[18], atan2(X[22], X[21]), sqrt(pow(X[21], 2) + pow(X[22], 2)));
			deltaAlf_pr1 = -promah * 15 * pi / 180 / 200;
			if (-deltaAlf_pr1 * 100 * sqrt(pow(X[21], 2) + pow(X[22], 2)) / g + cos(0) > 6)
			{
				deltaAlf_pr1 = -(6 - cos(0)) * g / sqrt(pow(X[21], 2) + pow(X[22], 2)) / 100;
			}

			if (-deltaAlf_pr1 * 100 * sqrt(pow(X[21], 2) + pow(X[22], 2)) / g + cos(0) < -6)
			{
				deltaAlf_pr1 = -(-6 - cos(0)) * g / sqrt(pow(X[21], 2) + pow(X[22], 2)) / 100;
			}
		}
		double deltaAlf_pr2 = 0;
		if (r_vizV(X[4], X[5], X[6], X[24], X[25], X[26]) < 15500)
		{
			Miss(X[5], X[4], Q, V, X[25], X[24], atan2(X[28], X[27]), sqrt(pow(X[27], 2) + pow(X[28], 2)));
			deltaAlf_pr2 = -promah * 15 * pi / 180 / 200;
			if (-deltaAlf_pr2 * 100 * sqrt(pow(X[27], 2) + pow(X[28], 2)) / g + cos(0) > 6)
			{
				deltaAlf_pr2 = -(6 - cos(0)) * g / sqrt(pow(X[27], 2) + pow(X[28], 2)) / 100;
			}
			if (-deltaAlf_pr2 * 100 * sqrt(pow(X[27], 2) + pow(X[28], 2)) / g + cos(0) < -6)
			{
				deltaAlf_pr2 = -(-6 - cos(0)) * g / sqrt(pow(X[27], 2) + pow(X[28], 2)) / 100;
			}
		}

		double Vx_newTg1, Vy_newTg1, Vx_newTg2, Vy_newTg2;
		double Vx_Tg1 = X[21];
		double Vy_Tg1 = X[22];
		double Vx_Tg2 = X[27];
		double Vy_Tg2 = X[28];
		for (int i = 0; i < 6; i++)
		{
			M_sv_sk(deltaAlf_pr1, 0);
			Vx_newTg1 = Matrix[0][0] * Vx_Tg1 + Matrix[0][1] * Vy_Tg1;
			Vy_newTg1 = Matrix[1][0] * Vx_Tg1 + Matrix[1][1] * Vy_Tg1;
			Vx_Tg1 = Vx_newTg1;
			Vy_Tg1 = Vy_newTg1;
			M_sv_sk(deltaAlf_pr2, 0);
			Vx_newTg2 = Matrix[0][0] * Vx_Tg2 + Matrix[0][1] * Vy_Tg2;
			Vy_newTg2 = Matrix[1][0] * Vx_Tg2 + Matrix[1][1] * Vy_Tg2;
			Vx_Tg2 = Vx_newTg2;
			Vy_Tg2 = Vy_newTg2;
			y_new_tg1 += Vy_Tg1;
			x_new_tg1 += Vx_Tg1;
			y_new_tg2 += Vy_Tg2;
			x_new_tg2 += Vx_Tg2;
		}
		x_tg = closer(x_new_tg1, x_new_tg2, x_new_tg1, x_new_tg2);
		y_tg = (y_new_tg1 + y_new_tg2) / 2;
		mu_new = atan2(y_tg - X[5], x_tg - X[4]) - Q;
	//////////////////////////////
	

/*	if (R_viz_12 >= 8000)
	{
		Xc = x_tg;
		Yc = y_tg;
		Zc = 0;
		Vxc = closer(x_new_tg1, x_new_tg2, Vx_newTg1, Vx_newTg2);
		Vyc = (Vy_newTg1 + Vy_newTg2) / 2;
		Vzc = 0;
		F[7] = r_vizV(X[4], X[5], X[6], Xc, Yc, Zc);
	}*/

/*
	V_Tg[0] = sqrt(pow(X[21], 2) + pow(X[22], 2));
	Theta_Tg[0] = atan2(X[22], X[21]);
	V_Tg[1] = sqrt(pow(X[27], 2) + pow(X[28], 2));
	Theta_Tg[1] = atan2(X[28], X[27]);
	V_Tg[2] = sqrt(pow(Vxc, 2) + pow(Vyc, 2));
	Vxc_12 = closer(X[18], X[24], V_Tg[0], V_Tg[1]);
	Theta_Tg[2] = theta_sr(closer(X[18], X[24], Theta_Tg[0], Theta_Tg[1]), closer(X[18], X[24], V_Tg[0], V_Tg[1]), closer(X[18], X[24], X[18], X[24]), further(X[18], X[24], Theta_Tg[0], Theta_Tg[1]), further(X[18], X[24], V_Tg[0], V_Tg[1]), further(X[18], X[24], X[19], X[25]));
	if (R_viz_12 >= 8000) //(target == 0)
						  //		else if (R_viz_12 >= 10000)
						  //		else if (t <= t_border - t_zad) //на втором участке траектории
	{
		Miss(X[5], X[4], Q, V, Yc, Xc, Theta_Tg[2], V_Tg[2]);
	}
	else //на участке траектории с конкретной целью
		 //			if (R_viz_12 >= 30) 
	{
		if (target == 1)
		{
			Miss(X[5], X[4], Q, V, Yc, Xc, Theta_Tg[0], V_Tg[0]);
		}
		if (target == 2)
		{
			Miss(X[5], X[4], Q, V, Yc, Xc, Theta_Tg[1], V_Tg[1]);
		}
	}
*/
/*	if (R_viz_12 >= 8000)
	{
	Xc = x_tg;
	Yc = y_tg;
	Zc = 0;
	Vxc = closer(x_new_tg1, x_new_tg2, Vx_newTg1, Vx_newTg2);
	Vyc = (Vy_newTg1 + Vy_newTg2) / 2;
	Vzc = 0;
	F[7] = r_vizV(X[4], X[5], X[6], Xc, Yc, Zc);
	}*/

	//углы наклона лиии визирования
	double fi_2 = viz_fi(X[5], F[7], Yc);
	double hi_2 = viz_hi(X[4], Xc, X[6], Zc);
//	if (Zc - X[6] > 0) K_hi *= -1;
	//матрица перевода из стартовой СК в сферическую
//	B_sf_st(fi_2, hi_2);
	B[0][0] = cos(fi_2)*cos(hi_2);
	B[0][1] = sin(fi_2);
	B[0][2] = -cos(fi_2)*sin(hi_2);
	B[1][0] = -sin(fi_2)*cos(hi_2);
	B[1][1] = cos(fi_2);
	B[1][2] = sin(fi_2)*sin(hi_2);
	B[2][0] = -sin(hi_2);
	B[2][1] = 0;
	B[2][2] = -cos(hi_2);

	//перевод взаимных скоростей цели и ракеты из стартовой СК в сферическую
	double Vr	= B[0][0] * (Vxc - Vxg) + B[0][1] * (Vyc - Vyg) + B[0][2] * (Vzc - Vzg);
	double Vfi	= B[1][0] * (Vxc - Vxg) + B[1][1] * (Vyc - Vyg) + B[1][2] * (Vzc - Vzg);
	double Vhi	= B[2][0] * (Vxc - Vxg) + B[2][1] * (Vyc - Vyg) + B[2][2] * (Vzc - Vzg);
//	Vhi = -sin(hi_2)*cos(fi_2) * (Vxc - Vxg) + cos(hi_2) * (Vyc - Vyg) + sin(hi_2)*sin(fi_2) * (Vzc - Vzg);
	//производные углов наклона линии визирования
	double d_fi = Vfi / F[7];
	double d_hi = Vhi / (F[7] * cos(fi_2));
//	d_hi = -Vhi / F[7];
	F[0] = fi_2;
	F[1] = hi_2;

	//производные углов наклона траектории ракеты
	double d_tet	= X[8] * sin(kren) + X[9] * cos(kren);
	double d_psi	= (1 / cos(tangaj))*(X[8] * cos(kren) - X[9] * sin(kren));
	double d_gamma	= X[8] - tan(tangaj)*(X[8] * cos(kren) - X[9] * sin(kren));

	//параметры отклонения рулей ракеты расчитанные в данный момент времени
	double del_v	= K1_v*K_fi*Kk_v*d_fi - K2_v*d_tet;
	double del_n	= K1_n*K_hi*Kk_n*d_hi - K2_n*d_psi;
//	double del_n	= +0.3*d_hi - 0.1*d_psi;
	double del_el	= - K1_el*d_gamma - K2_el*kren;

	double K3_v = -30;
	int sign = (epsilon > 0) - (epsilon < 0);
	double K4_v = sign * 45;

//	cout << sign << "\t" << epsilon * 180 / 3.141592 << "\n";

/*	if (R_viz_12 >= 8000)
	{
		del_v = K1_v*K_fi*Kk_v*d_fi - K2_v*d_tet + K3_v*epsilon + K4_v / F[7];
//		del_v = K1_v*K_fi*Kk_v*d_fi - K2_v*d_tet + K3_v*mu_new + K4_v / F[7];
	}*/

	

	//накладываемые ограничения на отклонения улей
	if (del_v > pi / 12)
		del_v = pi / 12;
	if (del_v < -pi / 12)
		del_v = -pi / 12;

	if (del_n > pi / 12)
		del_n = pi / 12;
	if (del_n < -pi / 12)
		del_n = -pi / 12;

	if (del_el > pi / 12)
		del_el = pi / 12;
	if (del_el < -pi / 12)
		del_el = -pi / 12;


//	if (tau >= t_gambit && (tau < t_gambit + 2) && exper == 0)
/*	if (tau >= 2 && (tau < 7) && exper == 0)
	{
		del_v = +pi / 20;
	}
	if (tau >= 7 && (tau < 10) && exper == 0)
	{
		del_v = -pi / 30;
	}
*/

	F[4] = del_v;
	F[5] = del_n;
	F[6] = del_el;

	//АД силы
	double FX = -Cx*q*S + P;
	double FY = Cy_a*q*S*alfa + Cydelta*q*S*del_v;
	double FZ = Cz_b*q*S*betta + Czdelta*q*S*del_n;
	//матрица переовода из стартовой СК в связанную
	RG_sv_st(X[10], X[11], X[12], X[13]);
	//переовод АД сил из связанной СК в стартовую
	double Xg = A[0][0] * FX + A[1][0] * FY + A[2][0] * FZ;
	double Yg = A[0][1] * FX + A[1][1] * FY + A[2][1] * FZ;
	double Zg = A[0][2] * FX + A[1][2] * FY + A[2][2] * FZ;

	double Ya = FX*sin(alfa) + FY*cos(alfa);
	double n_Ya = (P * sin(alfa) + Ya) / (m*g);	//нормальная перегрузка
	F[8] = n_Ya;
	//АД моменты
	double mx = mx_w*l*X[7] / V + Mxdel*del_el;
	double my = my_betta*betta + my_w*l*X[8] / V + Mydel*del_n;
	double mz = mz_alfa*alfa + mz_w*l*X[9] / V + Mzdel*del_v;

	//расчет правых частей парамтеров Родриго-Гамильтона
	IntRG(X[10], X[11], X[12], X[13], X[7], X[8], X[9]);
	//расчет правых частей ДУ
	Y0[1] = Xg / m; 		                                //Vxg
	Y0[2] = (Yg + G) / m;	                             	//Vyg
	Y0[3] = Zg / m;		                            		//Vzg
	Y0[4] = X[1];		                              		//Xg
	Y0[5] = X[2];		                              		//Yg
	Y0[6] = X[3];		                                	//Zg
	Y0[7] = mx / Ix + (Iy - Iz)*X[8] * X[9] / Ix;					//Omega x
	Y0[8] = my / Iy + (Iz - Ix)*X[7] * X[9] / Iy;					//Omega y
	Y0[9] = mz / Iz + (Ix - Iy)*X[7] * X[8] / Iz;					//Omega z

	Y0[10] = dRo;
	Y0[11] = dLambda;
	Y0[12] = dMu;
	Y0[13] = dNu;

	{
		deltaAlf_Tg[0] = 0;
	/*	if (r_vizV(X[4], X[5], X[6], X[18], X[19], X[20]) < 15000)
		{
			Miss(X[5], X[4], Q, V, X[19], X[18], atan2(X[22], X[21]), sqrt(pow(X[21], 2) + pow(X[22], 2)));
			deltaAlf_Tg[0] = (promah - 1200) * 15 * pi / 180 / 1500;
			if (-deltaAlf_Tg[0] * 100 * sqrt(pow(X[21], 2) + pow(X[22], 2)) / g + cos(-(atan2(X[22], X[21]) + pi)) > 5.6)
			{
				deltaAlf_Tg[0] = -(5.6 - cos(-(atan2(X[22], X[21]) + pi))) * g / sqrt(pow(X[21], 2) + pow(X[22], 2)) / 100;
			}

			if (-deltaAlf_Tg[0] * 100 * sqrt(pow(X[21], 2) + pow(X[22], 2)) / g + cos(-(atan2(X[22], X[21]) + pi)) < -4.0)
			{
				deltaAlf_Tg[0] = -(-4.0 - cos(-(atan2(X[22], X[21]) + pi))) * g / sqrt(pow(X[21], 2) + pow(X[22], 2)) / 100;
			}
			
		}*/
//		cout << (promah - 1200) * 15 / 1500 << "\t";
		if (X[19] > 10000 || X[19] < 2000)
		{
			deltaAlf_Tg[0] = -deltaAlf_Tg[0];
		}
		if (atan2(X[22], X[21]) > -160 * pi / 180 && atan2(X[22], X[21]) < 160 * pi / 180)
		{
			deltaAlf_Tg[0] = 0;
		}
		//матрица перевода из скоростной СК в связанную
		M_sv_sk(deltaAlf_Tg[0], deltaBet_Tg[0]);
		//изменение наклона траектории цели
		Vx_newTg[0] = Matrix[0][0] * X[21] + Matrix[0][1] * X[22] + Matrix[0][2] * X[23];
		Vy_newTg[0] = Matrix[1][0] * X[21] + Matrix[1][1] * X[22] + Matrix[1][2] * X[23];
		Vz_newTg[0] = Matrix[2][0] * X[21] + Matrix[2][1] * X[22] + Matrix[2][2] * X[23];
//		VxTG1 = Vx_newTg[0];
//		VyTG1 = Vy_newTg[0];
		VxTG1 = X[21];
		VyTG1 = X[22];

//		alfaTG = 15*(1 - tau / 12);

		deltaAlf_Tg[1] = 0;
	/*	if (r_vizV(X[4], X[5], X[6], X[24], X[25], X[26]) < 15000)
		{
			Miss(X[5], X[4], Q, V, X[25], X[24], atan2(X[28], X[27]), sqrt(pow(X[27], 2) + pow(X[28], 2)));
			deltaAlf_Tg[1] = (promah - 1200) * 15 * pi / 180 / 1500;
			if (-deltaAlf_Tg[1] * 100 * sqrt(pow(X[27], 2) + pow(X[28], 2)) / g + cos(-(atan2(X[28], X[27]) + pi)) > 5.6)
			{
				deltaAlf_Tg[1] = -(5.6 - cos(-(atan2(X[28], X[27]) + pi))) * g / sqrt(pow(X[27], 2) + pow(X[28], 2)) / 100;
			}
			if (-deltaAlf_Tg[1] * 100 * sqrt(pow(X[27], 2) + pow(X[28], 2)) / g + cos(-(atan2(X[28], X[27]) + pi)) < -4.0)
			{
				deltaAlf_Tg[1] = -(-4.0 - cos(-(atan2(X[28], X[27]) + pi))) * g / sqrt(pow(X[27], 2) + pow(X[28], 2)) / 100;
			}
		}*/
//		cout << (promah - 1200) * 15 / 1500 << "\n";
		if (X[25] > 10000 || X[25] < 2000)
		{
			deltaAlf_Tg[1] = -deltaAlf_Tg[1];
		}
		if (atan2(X[28], X[27]) > -160 * pi / 180 && atan2(X[28], X[27]) < 160 * pi / 180)
		{
			deltaAlf_Tg[1] = 0;
		}
		//матрица перевода из скоростной СК в связанную
		M_sv_sk(deltaAlf_Tg[1], deltaBet_Tg[1]);
		//изменение наклона траектории цели
		Vx_newTg[1] = Matrix[0][0] * X[27] + Matrix[0][1] * X[28] + Matrix[0][2] * X[29];
		Vy_newTg[1] = Matrix[1][0] * X[27] + Matrix[1][1] * X[28] + Matrix[1][2] * X[29];
		Vz_newTg[1] = Matrix[2][0] * X[27] + Matrix[2][1] * X[28] + Matrix[2][2] * X[29];
//		VxTG2 = Vx_newTg[1];
//		VyTG2 = Vy_newTg[1];
		VxTG2 = X[27];
		VyTG2 = X[28];
	}


//	cout << atan2(X[22], X[21]) * 180 / pi << "\t" << atan2(X[28], X[27]) * 180 / pi << "\n";
//	cout << -deltaAlf_Tg[0] * 100 * sqrt(pow(X[21], 2) + pow(X[22], 2)) / g + cos(-(atan2(X[22], X[21]) + pi)) << "\t" << -deltaAlf_Tg[1] * 100 * sqrt(pow(X[27], 2) + pow(X[28], 2)) / g + cos(-(atan2(X[28], X[27]) + pi)) << "\n";

	Y0[18] = Vx_newTg[0];
	Y0[19] = Vy_newTg[0];
	Y0[20] = Vz_newTg[0];
	Y0[24] = Vx_newTg[1];
	Y0[25] = Vy_newTg[1];
	Y0[26] = Vz_newTg[1];
/*
	Y0[21] = -10 +0.3 * tau;
	Y0[27] = -15;
*/	
//	cout << Xc << "\t"; cout << Yc << "\t"; cout << t << "\t"; cout << M << "\n";
	
	/*
	Y0[14] = d_gamma;
	Y0[15] = d_psi;
	Y0[16] = d_tet;
	*/
	if (R_viz_12 >= 8000) //(target == 0)
						  //		else if (R_viz_12 >= 10000)
						  //		else if (t <= t_border - t_zad) //на втором участке траектории
	{
		Miss(X[5], X[4], Q, V, Yc, Xc, Theta_Tg[2], V_Tg[2]);
	}
	else //на участке траектории с конкретной целью
		 //			if (R_viz_12 >= 30) 
	{
		if (target == 1)
		{
			Miss(X[5], X[4], Q, V, Yc, Xc, Theta_Tg[0], V_Tg[0]);
		}
		if (target == 2)
		{
			Miss(X[5], X[4], Q, V, Yc, Xc, Theta_Tg[1], V_Tg[1]);
		}
	}
}
//матрица переовода из стартовой СК в связанную через параметры Р-Г
void RG_sv_st(double Ro, double Lambda, double Mu, double Nu)
{
	A[0][0] = Ro*Ro + Lambda*Lambda - Mu*Mu - Nu*Nu;
	A[0][1] = 2 * (Ro*Nu + Lambda*Mu);
	A[0][2] = 2 * (-Ro*Mu + Lambda*Nu);

	A[1][0] = 2 * (-Ro*Nu + Lambda*Mu);
	A[1][1] = Ro*Ro + Mu*Mu - Nu*Nu - Lambda*Lambda;
	A[1][2] = 2 * (Ro*Lambda + Nu*Mu);

	A[2][0] = 2 * (Ro*Mu + Nu*Lambda);
	A[2][1] = 2 * (-Ro*Lambda + Nu*Mu);
	A[2][2] = Ro*Ro + Nu*Nu - Lambda*Lambda - Mu*Mu;
}
//матрица переовода из стартовой СК в связанную через углы Эйлера
void A_sv_st(double t, double r, double k)
{
	A[0][0] = cos(t)*cos(r);
	A[0][1] = sin(t);
	A[0][2] = -cos(t)*sin(r);

	A[1][0] = -sin(t)*cos(r)*cos(k) + sin(r)*sin(k);
	A[1][1] = cos(t)*cos(k);
	A[1][2] = sin(t)*sin(r)*cos(k) + cos(r)*sin(k);

	A[2][0] = sin(t)*cos(r)*sin(k) + sin(r)*cos(k);
	A[2][1] = -cos(t)*sin(k);
	A[2][2] = -sin(t)*sin(r)*sin(k) + cos(r)*cos(k);
}
//матрица перевода из стартовой СК в сферическую
void B_sf_st(double fi, double hi)
{
	B[0][0] = cos(fi)*cos(hi);
	B[0][1] = sin(fi);
	B[0][2] = -cos(fi)*sin(hi);
	B[1][0] = -sin(fi)*cos(hi);
	B[1][1] = cos(fi);
	B[1][2] = sin(fi)*sin(hi);
	B[2][0] = sin(hi);
	B[2][1] = 0;
	B[2][2] = cos(hi);
}
//матрица перевода из скоростной СК в связанную
void M_sv_sk(double alfa, double betta)
{
	Matrix[0][0] = cos(alfa)*cos(betta);
	Matrix[0][1] = -sin(alfa)*cos(betta);
	Matrix[0][2] = sin(betta);
	Matrix[1][0] = sin(alfa);
	Matrix[1][1] = cos(alfa);
	Matrix[1][2] = 0;
	Matrix[2][0] = -cos(alfa)*sin(betta);
	Matrix[2][1] = sin(alfa)*sin(betta);
	Matrix[2][2] = cos(betta);
};
//функция расчета скорсоти звука в атмосфере
double sk_zv(double H)
{
	double SoundSpeed;
	if (H <= 0.0)
	{
		SoundSpeed = 340.0;
	}
	else if (H > 0.0 && H <= 11000.0)
	{
		SoundSpeed = sqrt(299.79e-9*H*H - 2.6135*H + 115.80e+3);
	}
	else if (H > 11000.0 && H <= 25000.0)
	{
		SoundSpeed = 295.07;
	}
	else if (H > 25000.0 && H <= 46000.0)
	{
		SoundSpeed = sqrt(-179.80e-9*H*H + 1.1098*H + 59.433e+3);
	}
	else
	{
		SoundSpeed = 331.82;
	}
	return SoundSpeed;
};
//функция расчета плотности атмосферы
double pl(double H)
{
	double DensutyAir;
	if (H <= 0.0)
	{
		DensutyAir = 1.225;
	}
	else if (H > 0.0 && H <= 13000.0)
	{
		DensutyAir = 3.1088e-9*H*H - 114.04e-6*H + 1.223;
	}
	else if (H > 13000.0 && H <= 26000.0)
	{
		DensutyAir = log(-298.95e-12*H*H + (12.780e+3) / H + 28.677e-6*H);
	}
	else if (H > 26000.0 && H <= 39000.0)
	{
		DensutyAir = exp(-7.9310*log(H) - 94590 / H + 80.891);
	}
	else if (H > 39000.0 && H <= 49000.0)
	{
		DensutyAir = exp(1.4632e-9*H*H + 264.69e-3*log(H) - 266.74e-6*H);
	}
	else
	{
		DensutyAir = 0.0;
	}
	return DensutyAir;
};
//Функции нахождения АДХ атмосферы (включая производные по углам альфа, бетта, дельта)
double cx(double M, double alfa_)
{
	double v, n, koef;
	int j = -1, i = -1;
	//Зифа
	double m[10] = { 0.100, 0.500, 0.800, 0.900, 1.100, 1.200, 1.500, 1.900, 3.000, 100.000 };
	double al[5] = { 0.000, 2.000, 3.100, 5.000, 9.000 };
	double c[5][10] = { 0.521, 0.677, 0.747, 0.759, 0.824, 0.819, 0.803, 0.714, 0.549, 0.453,
						0.634, 0.678, 0.748, 0.759, 0.824, 0.820, 0.818, 0.727, 0.558, 0.460,
						0.637, 0.683, 0.753, 0.765, 0.826, 0.821, 0.850, 0.756, 0.578, 0.475,
						0.645, 0.693, 0.765, 0.776, 0.828, 0.824, 0.881, 0.783, 0.597, 0.491,
						0.672, 0.731, 0.806, 0.818, 0.837, 0.833, 0.950, 0.844, 0.642, 0.527 };

	alfa_ = alfa_ * 180 / pi;
	if (alfa_ < 0)
	{
		alfa_ = -alfa_;
	}
	if (alfa_ >= 8.0)
	{
		alfa_ = 7.999;
	}
	do
	{
		j++;
	} while (!(M >= m[j] && M<m[j + 1]));
	do
	{
		i++;
	}
	while (!(alfa_ >= al[i] && alfa_<al[i + 1]));
	n = (c[i][j + 1] - c[i][j])*(M - m[j]) / (m[j + 1] - m[j]) + c[i][j];
	v = (c[i + 1][j + 1] - c[i + 1][j])*(M - m[j]) / (m[j + 1] - m[j]) + c[i + 1][j];
	koef = (v - n)*(alfa_ - al[i]) / (al[i + 1] - al[i]) + n;
	return koef;
};

double cy_a(double M, double alfa_)
{
	double v, n, koef;
	int j = -1, i = -1;
	double m[10] = { 0.500, 0.900, 1.200, 1.500, 1.800, 2.100, 2.300, 2.500, 3.000, 100.000 };
	double al[5] = { 0.000, 1.000, 2.000, 3.000, 8.000 };
	double c[5][10] = { 1.250, 1.310, 1.418, 1.335, 1.266, 1.200, 1.150, 1.103, 0.979, 0.756,
						1.262, 1.331, 1.412, 1.329, 1.261, 1.196, 1.148, 1.102, 0.975, 0.746,
						1.289, 1.358, 1.412, 1.329, 1.261, 1.196, 1.148, 1.102, 0.975, 0.746,
						1.315, 1.386, 1.412, 1.329, 1.261, 1.196, 1.148, 1.102, 0.975, 0.746,
						1.447, 1.524, 1.412, 1.329, 1.261, 1.196, 1.148, 1.102, 0.975, 0.746 };

	alfa_ = alfa_ * 180 / pi;

	if (alfa_<0)
		alfa_ = -alfa_;

	if (alfa_ >= 8.0)
		alfa_ = 7.999;
	do
	{
		j++;
	}
	while (!(M >= m[j] && M<m[j + 1]));
	do
	{
		i++;
	}
	while (!(alfa_ >= al[i] && alfa_<al[i + 1]));
	n = (c[i][j + 1] - c[i][j])*(M - m[j]) / (m[j + 1] - m[j]) + c[i][j];
	v = (c[i + 1][j + 1] - c[i + 1][j])*(M - m[j]) / (m[j + 1] - m[j]) + c[i + 1][j];
	koef = (v - n)*(alfa_ - al[i]) / (al[i + 1] - al[i]) + n;
	return koef;
};

double cz_b(double M, double alfa_)
{
	return -cy_a(M, alfa_);
};

double cy_d(double M, double alfa_)
{
	double v, n, koef;
	int j = -1, i = -1;
	double m[10] = { 0.500, 0.900, 1.200, 1.500, 1.800, 2.100, 2.300, 2.500, 3.000, 100.000 };
	double al[5] = { 0.000, 1.000, 2.000, 3.000, 8.000 };
	double c[5][10] = { 0.116, 0.117, 0.118, 0.119, 0.119, 0.119, 0.119, 0.119, 0.119, 0.119,
						0.118, 0.119, 0.118, 0.119, 0.119, 0.119, 0.119, 0.119, 0.119, 0.119,
						0.120, 0.121, 0.118, 0.119, 0.119, 0.119, 0.119, 0.119, 0.119, 0.119,
						0.122, 0.123, 0.118, 0.119, 0.119, 0.119, 0.119, 0.119, 0.119, 0.119,
						0.133, 0.134, 0.118, 0.119, 0.119, 0.119, 0.119, 0.119, 0.119, 0.119 };

	alfa_ = alfa_ * 180 / pi;

	if (alfa_ < 0)
	{
		alfa_ = -alfa_;
	}
	if (alfa_ >= 8.0)
	{
		alfa_ = 7.999;
	}
	do
	{
		j++;
	}
	while (!(M >= m[j] && M<m[j + 1]));
	do
	{
		i++;
	}
	while (!(alfa_ >= al[i] && alfa_<al[i + 1]));
	n = (c[i][j + 1] - c[i][j])*(M - m[j]) / (m[j + 1] - m[j]) + c[i][j];
	v = (c[i + 1][j + 1] - c[ i + 1][j])*(M - m[j]) / (m[j + 1] - m[j]) + c[i + 1][j];
	koef = (v - n)*(alfa_ - al[i]) / (al[i + 1] - al[i]) + n;
	return koef*1.1;
};

double cz_d(double M, double alfa_)
{
	return -cy_d(M, alfa_);
};
double mx_wx(double M, double alfa_)
{
	return -0.005*0.6786;
};
double mz_wz(double M, double alfa_)
{
	double v, n, koef;
	int j = -1, i = -1, zn = 1;
	double m[10] = { 0.500, 0.900, 1.200, 1.500, 1.800, 2.100, 2.300, 2.500, 3.000, 100.000 };
	double al[5] = { 0.000, 1.000, 2.000, 3.000, 8.000 };
	double c[5][10] = { -3.281, -3.035, -3.480, -3.789, -3.727, -3.525, -3.398, -3.281, -2.990, -2.399,
						-3.357, -3.106, -3.480, -3.789, -3.727, -3.525, -3.398, -3.281, -2.990, -2.399,
						-3.432, -3.177, -3.480, -3.789, -3.727, -3.525, -3.398, -3.281, -2.990, -2.399,
						-3.508, -3.249, -3.480, -3.789, -3.727, -3.525, -3.398, -3.281, -2.990, -2.399,
						-3.888, -3.605, -3.480, -3.789, -3.727, -3.525, -3.398, -3.281, -2.990, -2.399 };
	alfa_ = alfa_ * 180 / pi;

	if (alfa_ < 0)
	{
		alfa_ = -alfa_;
	}
	if (alfa_ >= 8.0)
	{
		alfa_ = 7.999;
	}
	do
	{
		j++;
	}
	while (!(M >= m[j] && M<m[j + 1]));
	do
	{
		i++;
	}
	while (!(alfa_ >= al[i] && alfa_<al[i + 1]));
	n = (c[i][j + 1] - c[i][j])*(M - m[j]) / (m[j + 1] - m[j]) + c[i][j];
	v = (c[i + 1][j + 1] - c[i + 1][j])*(M - m[j]) / (m[j + 1] - m[j]) + c[i + 1][j];
	koef = (v - n)*(alfa_ - al[i]) / (al[i + 1] - al[i]) + n;
	return koef;
};
double mz_a(double M, double alfa_)
{
	double v, n, koef;
	int j = -1, i = -1;
	double m[10] = { 0.500, 0.900, 1.200, 1.500, 1.800, 2.100, 2.300, 2.500, 3.000, 100.000 };
	double al[5] = { 0.000, 1.000, 2.000, 3.000, 8.000 };
	double c[5][10] = { -0.099, -0.084, -0.095, -0.104, -0.101, -0.096, -0.092, -0.089, -0.078, -0.054,
						-0.098, -0.084, -0.094, -0.103, -0.100, -0.095, -0.092, -0.089, -0.077, -0.052,
						-0.099, -0.085, -0.094, -0.103, -0.100, -0.095, -0.092, -0.089, -0.077, -0.052,
						-0.101, -0.086, -0.094, -0.103, -0.100, -0.095, -0.092, -0.089, -0.077, -0.052,
						-0.108, -0.091, -0.094, -0.103, -0.100, -0.095, -0.092, -0.089, -0.077, -0.052 };

	alfa_ = alfa_ * 180 / pi;

	if (alfa_ < 0)
	{
		alfa_ = -alfa_;
	}
	if (alfa_ >= 8.0)
	{
		alfa_ = 7.999;
	}
	do
	{
		j++;
	}
	while (!(M >= m[j] && M<m[j + 1]));
	do
	{
		i++;
	}
	while (!(alfa_ >= al[i] && alfa_<al[i + 1]));
	n = (c[i][j + 1] - c[i][j])*(M - m[j]) / (m[j + 1] - m[j]) + c[i][j];
	v = (c[i + 1][j + 1] - c[i + 1][j])*(M - m[j]) / (m[j + 1] - m[j]) + c[i + 1][j];
	koef = (v - n)*(alfa_ - al[i]) / (al[i + 1] - al[i]) + n;
	return koef;
};

double mz_d(double M, double alfa_)
{
	double v, n, koef;
	int j = -1, i = -1;
	double m[10] = { 0.500, 0.900, 1.200, 1.500, 1.800, 2.100, 2.300, 2.500, 3.000, 100.000 };
	double al[5] = { 0.000, 1.000, 2.000, 3.000, 8.000 };
	double c[5][10] = { 0.081, 0.083, 0.090, 0.093, 0.095, 0.097, 0.097, 0.097, 0.097, 0.095,
						0.091, 0.094, 0.090, 0.093, 0.095, 0.097, 0.097, 0.097, 0.097, 0.095,
						0.102, 0.104, 0.090, 0.093, 0.095, 0.097, 0.097, 0.097, 0.097, 0.095,
						0.112, 0.115, 0.090, 0.093, 0.095, 0.097, 0.097, 0.097, 0.097, 0.095,
						0.163, 0.168, 0.090, 0.093, 0.095, 0.097, 0.097, 0.097, 0.097, 0.095 };

	alfa_ = alfa_ * 180 / pi;

	if (alfa_ < 0)
	{
		alfa_ = -alfa_;
	}
	if (alfa_ >= 8.0)
	{
		alfa_ = 7.999;
	}
	do
	{
		j++;
	}
	while (!(M >= m[j] && M<m[j + 1]));
	do
	{
		i++;
	}
	while (!(alfa_ >= al[i] && alfa_<al[i + 1]));
	n = (c[i][j + 1] - c[i][j])*(M - m[j]) / (m[j + 1] - m[j]) + c[i][j];
	v = (c[i + 1][j + 1] - c[i + 1][j])*(M - m[j]) / (m[j + 1] - m[j]) + c[i + 1][j];
	koef = (v - n)*(alfa_ - al[i]) / (al[i + 1] - al[i]) + n;
	return koef;
}; 
double my_wy(double M, double alfa_)
{
	return mz_wz(M, alfa_);
};
double my_b(double M, double alfa_)
{
	return mz_a(M, alfa_);
};
double my_d(double M, double alfa_)
{
	return mz_d(M, alfa_);
};
double mx_d(double M, double alfa_)
{
	return 0.1*mz_d(M, alfa_);
};
//Функция определния наклона линии визирования в вертикальной плоскости
double viz_fi(double Y, double r_viz, double Yc)
{
	return asin((Yc - Y) / r_viz);
};

//Функция определния наклона линии визирования в горизонтальной плоскости
double viz_hi(double X, double Xc, double Z, double Zc)
{
	return -atan2(Zc - Z, Xc - X);
//	return -atan((Zc - Z) / (Xc - X));

};

//Функция определния радиуса линии визирования
double r_vizV(double X, double Y, double Z, double Xc, double Yc, double Zc)
{
	return sqrt(pow((Yc - Y), 2) + pow((Xc - X), 2) + pow((Zc - Z), 2));
};

void Miss(double Y_LA, double X_LA, double theta_la, double V_la, double Y_TG, double X_TG, double theta_tg, double V_tg)
{
	double r_v = sqrt(pow((Y_TG - Y_LA), 2) + pow((X_TG - X_LA), 2));
	double y_la = Y_LA, x_la = X_LA;
	y_tg = Y_TG;
	x_tg = X_TG;
	do
	{
		promah = r_v;
		y_la += 0.01 * V_la * sin(theta_la);
		x_la += 0.01 * V_la * cos(theta_la);
		y_tg += 0.01 * V_tg * sin(theta_tg);
		x_tg += 0.01 * V_tg * cos(theta_tg);
		r_v = sqrt(pow((y_tg - y_la), 2) + pow((x_tg - x_la), 2));
	} while (promah > r_v);
	y_la -= 0.01 * V_la * sin(theta_la);
	x_la -= 0.01 * V_la * cos(theta_la);
	y_tg -= 0.01 * V_tg * sin(theta_tg);
	x_tg -= 0.01 * V_tg * cos(theta_tg);
	epsilon = atan2(y_tg - Y_LA, x_tg - X_LA) - atan2(y_la - Y_LA, x_la - X_LA);
	promah *= (epsilon > 0) - (epsilon < 0);
	promah *= -1;
};

double theta_sr(double theta_closer, double V_closer, double xn_closer, double theta_further, double V_further, double yn_further)
{
	double xk_closer = xn_closer + V_closer * cos(theta_closer);
	double yk_further = yn_further + V_further * cos(theta_further);
	return (theta_closer + atan2(yk_further - yn_further, xk_closer - xn_closer)) / 2;
}

void Prognoz(double y_new_tg, double x_new_tg, double Vy_new_tg, double Vx_new_tg, double Y_LA, double X_LA, double theta_la)
{
	y_new_tg += 5 * Vy_new_tg;
	x_new_tg += 5 * Vx_new_tg;
	mu_new = atan2(y_new_tg - Y_LA, x_new_tg - X_LA) - theta_la;
}
/*
template <typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}
*/