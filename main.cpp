#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
//Методы оптимизации: одномерная минимизация
//Функция: e^(6*x) + e^(-7*x) + 7*x^2 + 0*x + 6 - x^3

std::ofstream fout("output.txt");

struct function {
    static int counter, derCounter, secderCounter;

    static double func(double x) {
        //Возвращает значение функции f(x) в точке x
        counter++; //Счетчик вычислений функций
        return exp(6 * x) + exp(-7 * x) + 7 * pow(x, 2) + 6 - pow(x, 3);
    }

    static double der(double x) {
        //Возвращает значение производной функции f'(x) в точке x
        derCounter++; //Счетчик вычислений функций
        return 6 * exp(6 * x) - 7 * exp(-7 * x) + 14 * x - 3 * pow(x, 2);
    }

    static double secDer(double x) {
        //Возвращает значение второй производной функции f"(x) в точке x
        secderCounter++; //Счетчик вычислений функций
        return 36 * exp(6 * x) + 49 * exp(-7 * x) + 14 - 6 * x;
    }

    static void printResult(double x, double y){
        fout << "• Минимум в точке: " << x << ", и он равен: " << y /*<< std::setprecision(6)*/ << std::endl;
        fout << "• Приближенные минимайзеры: " << std::setprecision(4) << x << ", " << y;
        fout << " (функция была вычислена " << counter << " раз(а), её производная — " << derCounter << " раз(а), "
                                        "вторая производная — " << secderCounter << " раз(а))"<< std::endl;
    }

    static void printError(double xleft, double xright, double yleft, double yright) {
        //Вычисляет абсолютную погрешность
        fout << "• Δx = " << abs(xright - xleft) / 2 << ", Δy = " << abs(yright - yleft) / 2 << std::endl;
    }
};
int function::counter = 0, function::derCounter = 0, function::secderCounter = 0;


void PassiveSearch(double epsilon, double left, double right) {
    //Метод пассивного поиска
    fout << "МЕТОД ПАССИВНОГО ПОИСКА" << std::endl;
    fout << std::setprecision(15);
    int k = ceil((right - left) / epsilon);
    double xmin = left, ymin = function::func(left);
    fout << "┌────────┬─────────┬──────────────────┐\n"
            "│   i    │    x    │       f(x)       │\n"
            "├────────┼─────────┼──────────────────┤\n"
            ;
    fout << std::left << "│ 0      " << "│ " << std::setw(8) << xmin << "│ " << std::setw(17) << ymin << "│" << std::endl;
    for(int i = 0; i < k; i++) { //Критерий остановки
        double xi = left + (right - left) * (i + 1) / k;
        double yi = function::func(xi);
        fout << "│ " << std::setw(7) << i + 1 << "│ " << std::setw(8) << xi;
        fout << "│ " << std::setw(17) << yi << "│" << std::endl;
        if (yi < ymin) {
            xmin = xi; ymin = yi;
        } else
            break; //т.к. наша функция унимодальна
    }
    fout << "└────────┴─────────┴──────────────────┘" << std::endl;
    fout << "• Здесь x(min) = min(x(i))" << std::endl;
    function::printResult(xmin, ymin);
    function::counter = 0; //Обнуляем счетчик
    fout << std::endl;
}

void Dichotomy(double epsilon, double left, double right) {
    //Метод дихотомии
    fout << "МЕТОД ДИХОТОМИИ" << std::endl;
    fout << std::setprecision(15);
    int i = 0;
    double delta = epsilon / 2;
    double ai = left, bi = right, ci, di;
    double FuncCi, FuncDi;
    fout << "┌────────┬──────────────────────┬──────────────────────┬──────────────────────┬──────────────────────┬──────────────────────┬──────────────────────┐\n"
            "│   i    │         a(i)         │         b(i)         │         c(i)         │         d(i)         │         F(c)         │         F(d)         │\n"
            "├────────┼──────────────────────┼──────────────────────┼──────────────────────┼──────────────────────┼──────────────────────┼──────────────────────┤\n"
            ;
    while((bi - ai) / 2 > epsilon) { //Критерий остановки
        ci = (ai + bi) / 2 - delta / 2;
        di = (ai + bi) / 2 + delta / 2;
        FuncCi = function::func(ci); FuncDi = function::func(di);
        fout << std::left;
        fout << "│ " << std::setw(7) << i << "│ " << std::setw(21)<< ai << "│ " << std::setw(21) << bi;
        fout << "│ " << std::setw(21) << ci << "│ " << std::setw(21)<< di;
        fout << "│ " << std::setw(21) << FuncCi << "│ " << std::setw(21)<< FuncDi << "│ " << std::endl;
        if(FuncCi <= FuncDi) bi = di;
        else ai = ci;
        i++;
    }
    fout << "│ " << std::setw(7) << i << "│ " << std::setw(21)<< ai << "│ " << std::setw(21) << bi;
    fout << "│ " << std::setw(21) << ci << "│ " << std::setw(21)<< di;
    fout << "│ " << std::setw(21) << FuncCi << "│ " << std::setw(21)<< FuncDi << "│ " << std::endl;
    fout << "└────────┴──────────────────────┴──────────────────────┴──────────────────────┴──────────────────────┴──────────────────────┴──────────────────────┘";
    fout << std::endl;
    fout << "• Здесь x(min) = (a(i) + b(i)) / 2" << std::endl;
    double xmin = (ai + bi) / 2;
    function::printResult(xmin, function::func(xmin));
    function::printError(ai, bi, function::func(ai), function::func(bi));
    function::counter = 0; //Обнуляем счетчик
    fout << std::endl;
}

void GoldenRatio(double epsilon, double left, double right){
    //Метод золотого сечения
    fout << "МЕТОД ЗОЛОТОГО СЕЧЕНИЯ" << std::endl;
    fout << std::setprecision(15);
    int i = 0;
    double ai = left, bi = right;
    fout << "┌────────┬──────────────────────┬──────────────────────┬──────────────────────┬──────────────────────┬──────────────────────┬──────────────────────┐\n"
            "│   i    │         a(i)         │         b(i)         │         c(i)         │         d(i)         │         F(c)         │         F(d)         │\n"
            "├────────┼──────────────────────┼──────────────────────┼──────────────────────┼──────────────────────┼──────────────────────┼──────────────────────┤\n"
            ;
    double ci = (3 - sqrt(5)) * (bi -ai) / 2 + ai;
    double di = (sqrt(5) - 1) * (bi -ai) / 2 + ai;
    double Fci = function::func(ci), Fdi = function::func(di);
    while((bi - ai) / 2 > epsilon) { //Критерий остановки
        fout << std::left;
        fout << "│ " << std::setw(7) << i << "│ " << std::setw(21)<< ai << "│ " << std::setw(21) << bi;
        fout << "│ " << std::setw(21) << ci << "│ " << std::setw(21)<< di;
        fout << "│ " << std::setw(21) << Fci << "│ " << std::setw(21)<< Fdi << "│ " << std::endl;
        if(Fci <= Fdi) {
            bi = di;
            di = ci;
            ci = (3 - sqrt(5)) * (bi -ai) / 2 + ai;
            Fdi = Fci;
            Fci = function::func(ci);
        } else {
            ai = ci;
            ci = di;
            di = (sqrt(5) - 1) * (bi -ai) / 2 + ai;
            Fci = Fdi;
            Fdi = function::func(di);
        }
        i++;
    }
    fout << "│ " << std::setw(7) << i << "│ " << std::setw(21)<< ai << "│ " << std::setw(21) << bi;
    fout << "│ " << std::setw(21) << ci << "│ " << std::setw(21)<< di;
    fout << "│ " << std::setw(21) << Fci << "│ " << std::setw(21)<< Fdi << "│ " << std::endl;
    fout << "└────────┴──────────────────────┴──────────────────────┴──────────────────────┴──────────────────────┴──────────────────────┴──────────────────────┘";
    fout << std::endl;
    fout << "• Здесь x(min) = (a(i) + b(i)) / 2" << std::endl;
    double xmin = (ai + bi) / 2;
    function::printResult(xmin, function::func(xmin));
    function::printError(ai, bi, function::func(ai), function::func(bi));
    function::counter = 0; //Обнуляем счетчик
    fout << std::endl;
}

void Fibonacci(double epsilon, double left, double right) {
    //Метод Фибоначчи (упрощенный)
    fout << "МЕТОД ФИБОНАЧЧИ" << std::endl;
    fout << std::setprecision(15);
    int t = 1;
    double ai = left, bi = right;
    std::vector<double> F;
    F.assign(2, 1);
    //Находим n, если учесть, что F(n+2) >= (b-a)/ε
    while((bi - ai) / epsilon > F[t]) {
        t++;
        F.push_back(F[t - 1] + F[t - 2]);
    }
    int n = t - 2;
    fout << "┌────────┬──────────────────────┬──────────────────────┬──────────────────────┬──────────────────────┬──────────────────────┬──────────────────────┐\n"
            "│   j    │         a(j)         │         b(j)         │         c(j)         │         d(j)         │         F(c)         │         F(d)         │\n"
            "├────────┼──────────────────────┼──────────────────────┼──────────────────────┼──────────────────────┼──────────────────────┼──────────────────────┤\n"
            ;
    double ci = ai + (bi - ai) * (F[n] / F[n + 2]);
    double di = ai + (bi - ai) * (F[n + 1] / F[n + 2]);
    double Fci = function::func(ci), Fdi = function::func(di);
    for(int j = 0; j < n; j++) { //Критерий остановки
        fout << std::left;
        fout << "│ " << std::setw(7) << j << "│ " << std::setw(21)<< ai << "│ " << std::setw(21) << bi;
        fout << "│ " << std::setw(21) << ci << "│ " << std::setw(21)<< di;
        fout << "│ " << std::setw(21) << Fci << "│ " << std::setw(21)<< Fdi << "│ " << std::endl;
        if(Fci <= Fdi) {
            bi = di;
            di = ci;
            ci = ai + (bi - ai) * (F[n - j] / F[n + 2 - j]);
            Fdi = Fci;
            Fci = function::func(ci);
        } else {
            ai = ci;
            ci = di;
            di = ai + (bi - ai) * (F[n - j + 1] / F[n + 2 - j]);
            Fci = Fdi;
            Fdi = function::func(di);
        }
    }
    fout << "│ " << std::setw(7) << n << "│ " << std::setw(21)<< ai << "│ " << std::setw(21) << bi;
    fout << "│ " << std::setw(21) << ci << "│ " << std::setw(21)<< di;
    fout << "│ " << std::setw(21) << Fci << "│ " << std::setw(21)<< Fdi << "│ " << std::endl;
    fout << "└────────┴──────────────────────┴──────────────────────┴──────────────────────┴──────────────────────┴──────────────────────┴──────────────────────┘";
    fout << std::endl;
    fout << "• Здесь x(min) = (a(j) + b(j)) / 2" << std::endl;
    double xmin = (ai + bi) / 2;
    function::printResult(xmin, function::func(xmin));
    function::printError(ai, bi, function::func(ai), function::func(bi));
    function::counter = 0; //Обнуляем счетчик
    fout << std::endl;
}

void Tangent(double epsilon, double left, double right) {
    //Метод касательных
    fout << "МЕТОД КАСАТЕЛЬНЫХ" << std::endl;
    fout << std::setprecision(15);
    int i = 0;
    double ai = left, bi = right, ci, xmin, Fxmin;
    double FuncAi = function::func(ai), FuncBi = function::func(bi), FuncCi;
    double DerAi = function::der(ai), DerBi = function::der(bi), DerCi;
    if (DerAi >= 0) {
        xmin = ai;
        Fxmin = FuncAi;
    } else if (DerBi <= 0) {
        xmin = bi;
        Fxmin = FuncBi;
    } else {
        fout << "┌────────┬──────────────────────┬──────────────────────┬──────────────────────┬──────────────────────┬───────────────────────┐\n"
                "│   i    │         a(i)         │         b(i)         │         c(i)         │         F(c)         │         F′(c)         │\n"
                "├────────┼──────────────────────┼──────────────────────┼──────────────────────┼──────────────────────┼───────────────────────┤\n"
                ;
        while (abs(bi - ai) > epsilon) {
            ci = (FuncBi - FuncAi + DerAi * ai - DerBi * bi) / (DerAi - DerBi);
            FuncCi = function::func(ci);
            DerCi = function::der(ci);
            fout << std::left;
            fout << "│ " << std::setw(7) << i << "│ " << std::setw(21)<< ai << "│ " << std::setw(21) << bi;
            fout << "│ " << std::setw(21) << ci;
            fout << "│ " << std::setw(21)<< FuncCi << "│ " << std::setw(22)<< DerCi << "│" << std::endl;
            if (DerCi < 0) {
                ai = ci;
                FuncAi = FuncCi;
                DerAi = DerCi;
            } else if (DerCi > 0) {
                bi = ci;
                FuncBi = FuncCi;
                DerBi = DerCi;
            } else break;
            i++;
        }
        xmin = ci;
        Fxmin = FuncCi;
        fout << "│ " << std::setw(7) << i << "│ " << std::setw(21)<< ai << "│ " << std::setw(21) << bi;
        fout << "│ " << std::setw(21) << ci;
        fout << "│ " << std::setw(21)<< FuncCi << "│ " << std::setw(22)<< DerCi << "│" << std::endl;
        fout << "└────────┴──────────────────────┴──────────────────────┴──────────────────────┴──────────────────────┴───────────────────────┘";
        fout << std::endl;
    }
    fout << "• Здесь x(min) = c(i)" << std::endl;
    function::printResult(xmin, Fxmin);
    function::counter = 0, function::derCounter = 0; //Обнуляем счетчики
    fout << std::endl;
}

void NewtonRaphson(double epsilon, double xi) {
    //Метод Ньютона-Рафсона
    fout << "МЕТОД НЬЮТОНА-РАФСОНА" << std::endl;
    fout << std::setprecision(15);
    int i = 0;
    fout << "┌────────┬──────────────────────┬──────────────────────┬──────────────────────┐\n"
            "│   i    │         x(i)         │         f′(x)        │         f″(x)        │\n"
            "├────────┼──────────────────────┼──────────────────────┼──────────────────────┤\n"
            ;
    double DerXi = function::der(xi), SecDerXi = function::secDer(xi);
    while (abs(DerXi) > epsilon) {
        fout << std::left;
        fout << "│ " << std::setw(7) << i << "│ " << std::setw(21)<< xi;
        fout << "│ " << std::setw(21) << DerXi << "│ " << std::setw(21)<< SecDerXi << "│" << std::endl;
        xi = xi - (DerXi / SecDerXi);
        DerXi = function::der(xi), SecDerXi = function::secDer(xi);
        i++;
    }
    fout << "│ " << std::setw(7) << i << "│ " << std::setw(21)<< xi;
    fout << "│ " << std::setw(21) << DerXi << "│ " << std::setw(21)<< SecDerXi << "│" << std::endl;
    fout << "└────────┴──────────────────────┴──────────────────────┴──────────────────────┘" << std::endl;
    fout << "• Здесь x(min) = x(i)" << std::endl;
    function::printResult(xi, function::func(xi));
    function::counter = 0, function::derCounter = 0, function::secderCounter = 0; //Обнуляем счетчики
    fout << std::endl;
}

//void Secant(double epsilon, double x1, double x2) {
//
//}


int main() {
//    PassiveSearch(0.001, 0, 0.102); //Здесь точность ниже, чем нужно, чтобы файл вывода не весил очень много
    PassiveSearch(0.00001, 0, 0.102);
    Dichotomy(0.00001, 0, 1);
    GoldenRatio(0.00001, 0, 1);
    Fibonacci(0.00001, 0, 1);
    Tangent(0.000001, 0, 1);
    NewtonRaphson(0.00001, 0.5);
    fout.close();
    return 0;
}
