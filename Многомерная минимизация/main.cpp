#include <iostream>
#include <fstream>
#include <complex>
#include <optional>
#include <vector>
//Методы оптимизации: многомерная минимизация
//Функция (по z): (z - 7) * (z - 7i) * (z - 7 - 6i) = z^3 - (14 + 13i) * z^2 + (7 + 140i) * z + (294 - 343i)
//Функция (по x, y): x^3 + 3i*x^2*y - (14 + 13i)x^2 - 3xy^2 + (26 - 28i)xy + (7 + 140i)x - iy^3 + (14 + 13i)y^2 - (140 - 7i)y + (294 - 343i)
//Минимайзер (по х, у): (x^2 - 14x + y^2 + 49) * (x^2 + y^2 - 14y + 49) * (x^2 - 14x + y^2 - 12y + 85);

std::ofstream fout("output.txt");

struct function {
    static int counter, gradientCounter;

    static std::complex<double> func(double x, double y) {
        //Возвращает значение функции P(x, y) в точке (x, y)
        std::complex<double> i(0.0, 1.0);
        counter++; //Счетчик вычислений функций
        return pow(x, 3) + 3.0 * i * pow(x, 2) * y - (14.0 + 13.0 * i) * pow(x, 2) - 3.0 * x * pow(y, 2) + (26.0 - 28.0 * i) * x * y + (7.0 + 140.0 * i) * x - pow(y, 3) * (0.0 + 1.0 * i) + (14.0 + 13.0 * i) * pow(y, 2) - (140.0 - 7.0 * i) * y + (294.0 - 343.0 * i);
    }

    static double minimizer(double x, double y,
                            std::optional<double> x1 = std::nullopt, std::optional<double> y1 = std::nullopt,
                            std::optional<double> x2 = std::nullopt, std::optional<double> y2 = std::nullopt) {
        //Возвращает значение функции |P(x, y)|^2 в точке (x, y)
        counter++; //Счетчик вычислений функций
        double func = (pow(x, 2) - 14 * x + pow(y, 2) + 49) * (pow(x, 2) + pow(y, 2) - 14 * y + 49) * (pow(x, 2) - 14 * x + pow(y, 2) - 12 * y + 85);
        if (x1 && y1) {
            double denom1 = pow(x, 2) + pow(y, 2) + pow(*x1, 2) + pow(*y1, 2) + 2 * (- x * *x1 - y * *y1);
            if (x2 && y2) {
                double denom2 = pow(x, 2) + pow(y, 2) + pow(*x2, 2) + pow(*y2, 2) + 2 * (- x * *x2 - y * *y2);
                return func / (denom1 * denom2);
            }
            return func / denom1;
        }
        return func;
    }

    static bool stepAdvance(double& x, double& y, double& func,
                            double step_x, double step_y, double true_step_x, double true_step_y,
                            int& i, const std::string& type,
                            std::optional<double> x1 = std::nullopt, std::optional<double> y1 = std::nullopt,
                            std::optional<double> x2 = std::nullopt, std::optional<double> y2 = std::nullopt) {
        //Минимизация и печать для покоординатного спуска
        fout << "│ " << std::setw(6) << type;
        fout << "│ " << std::setw(21)<< true_step_x << "│ " << std::setw(21) << true_step_y;
        fout << "│ " << std::setw(20)<< (i + 1) << "│ " << std::endl;
        i++;
        double new_x = x + step_x, new_y = y + step_y;
        double new_func = minimizer(new_x, new_y, x1, y1, x2, y2);
        if (new_func < func) {
            while (new_func < func) {
                x = new_x; y = new_y;
                func = new_func;
                fout << "│ " << std::setw(7) << i << "│ " << std::setw(21)<< x << "│ " << std::setw(21) << y;
                fout << "│ " << std::setw(21)<< func;
                fout << "│ " << std::setw(6) << type;
                fout << "│ " << std::setw(21)<< true_step_x << "│ " << std::setw(21) << true_step_y;
                fout << "│ " << std::setw(20)<< (i + 1) << "│ " << std::endl;
                i++;
                new_x += step_x; new_y += step_y;
                new_func = minimizer(new_x, new_y, x1, y1, x2, y2);
                if (new_func >= func) {
                    fout << "│ " << std::setw(7) << i << "│ " << std::setw(21) << new_x << "│ " << std::setw(21) << new_y;
                    fout << "│ " << std::setw(21)<< new_func;
                }
            }
            return true;
        } else {
            fout << "│ " << std::setw(7) << i << "│ " << std::setw(21)<< new_x << "│ " << std::setw(21) << new_y;
            fout << "│ " << std::setw(21)<< new_func;
            return false;
        }
    }

    static double alphaFunction(double x, double y, double alpha,
                                std::optional<double> x1 = std::nullopt, std::optional<double> y1 = std::nullopt,
                                std::optional<double> x2 = std::nullopt, std::optional<double> y2 = std::nullopt) {
        //Возвращает значение функции |P(x, y)|^2 в точке (x - αP′x, y - αP′y)
        counter++; //Счетчик вычислений функций
        double Der1 = function::xDerivative(x, y, x1, y1, x2, y2), Der2 = function::yDerivative(x, y, x1, y1, x2, y2);
        double func = (pow(x - alpha * Der1, 2) - 14 * (x - alpha * Der1) + pow((y - alpha * Der2), 2) + 49) * (pow(x - alpha * Der1, 2) + pow((y - alpha * Der2), 2) - 14 * (y - alpha * Der2) + 49) * (pow(x - alpha * Der1, 2) - 14 * (x - alpha * Der1) + pow((y - alpha * Der2), 2) - 12 * (y - alpha * Der2) + 85);
        if (x1 && y1) {
            double denom1 = pow(x - alpha * Der1, 2) + pow((y - alpha * Der2), 2) + pow(*x1, 2) + pow(*y1, 2) + 2 * (- (x - alpha * Der1) * *x1 - (y - alpha * Der2) * *y1);
            if (x2 && y2) {
                double denom2 = pow(x - alpha * Der1, 2) + pow((y - alpha * Der2), 2) + pow(*x2, 2) + pow(*y2, 2) + 2 * (- (x - alpha * Der1) * *x2 - (y - alpha * Der2) * *y2);
                return func / (denom1 * denom2);
            }
            return func / denom1;
        }
        return func;
    }

    static double xDerivative(double x, double y,
                              std::optional<double> x1 = std::nullopt, std::optional<double> y1 = std::nullopt,
                              std::optional<double> x2 = std::nullopt, std::optional<double> y2 = std::nullopt) {
        //Возвращает значение производной функции |P(x, y)|^2 по х в точке (x, y)
        gradientCounter++; //Счетчик вычислений производной
        double funcDer = 6*pow(x, 5) - 140*pow(x, 4) + 12*pow(x, 3)*pow(y, 2) - 104*pow(x, 3)*y + 1516*pow(x, 3) - 168*pow(x, 2)*pow(y, 2) + 1680*pow(x, 2)*y - 9744*pow(x, 2) + 6*x*pow(y, 4) - 104*x*pow(y, 3) + 1460*x*pow(y, 2) - 11592*x*y + 40670*x - 28*pow(y, 4) + 560*pow(y, 3) - 5600*pow(y, 2) + 34496*y - 91924;
        if (x1 && y1) {
            double func = (pow(x, 2) - 14 * x + pow(y, 2) + 49) * (pow(x, 2) + pow(y, 2) - 14 * y + 49) * (pow(x, 2) - 14 * x + pow(y, 2) - 12 * y + 85);
            double denom1 = pow(x, 2) + pow(y, 2) + pow(*x1, 2) + pow(*y1, 2) + 2 * (- x * *x1 - y * *y1);
            double denom1Der = 2 * (x - *x1);
            if (x2 && y2) {
                double denom2 = pow(x, 2) + pow(y, 2) + pow(*x2, 2) + pow(*y2, 2) + 2 * (- x * *x2 - y * *y2);
                double denom2Der = 2 * (x - *x2);
                return (funcDer * denom1 * denom2 - func * denom1Der * denom2 - func * denom1 * denom2Der) / (pow(denom1, 2) * pow(denom2, 2));
            }
            return (funcDer * denom1 - func * denom1Der) / pow(denom1, 2);
        }
        return funcDer;
    }

    static double yDerivative(double x, double y,
                              std::optional<double> x1 = std::nullopt, std::optional<double> y1 = std::nullopt,
                              std::optional<double> x2 = std::nullopt, std::optional<double> y2 = std::nullopt) {
        //Возвращает значение производной функции |P(x, y)|^2 по y в точке (x, y)
        gradientCounter++; //Счетчик вычислений функций
        double funcDer = 6*pow(x, 4)*y - 26*pow(x, 4) - 112*pow(x, 3)*y + 560*pow(x, 3) + 12*pow(x, 2)*pow(y, 3) - 156*pow(x, 2)*pow(y, 2) + 1460*pow(x, 2)*y - 5796*pow(x, 2) - 112*x*pow(y, 3) + 1680*x*pow(y, 2) - 11200*x*y + 34496*x + 6*pow(y, 5) - 130*pow(y, 4) + 1404*pow(y, 3) - 9156*pow(y, 2) + 37926*y - 87122;
        if (x1 && y1) {
            double func = (pow(x, 2) - 14 * x + pow(y, 2) + 49) * (pow(x, 2) + pow(y, 2) - 14 * y + 49) * (pow(x, 2) - 14 * x + pow(y, 2) - 12 * y + 85);
            double denom1 = pow(x, 2) + pow(y, 2) + pow(*x1, 2) + pow(*y1, 2) + 2 * (- x * *x1 - y * *y1);
            double denom1Der = 2 * (y - *y1);
            if (x2 && y2) {
                double denom2 = pow(x, 2) + pow(y, 2) + pow(*x2, 2) + pow(*y2, 2) + 2 * (- x * *x2 - y * *y2);
                double denom2Der = 2 * (y - *y2);
                return (funcDer * denom1 * denom2 - func * denom1Der * denom2 - func * denom1 * denom2Der) / (pow(denom1, 2) * pow(denom2, 2));
            }
            return (funcDer * denom1 - func * denom1Der) / pow(denom1, 2);
        }
        return funcDer;
    }

    static void printResult(double x, double y, double z){
        fout << "• Минимум в точке: (" << x << ", " << y << "), и он равен: " << z << std::endl;
        fout << "• Приближенные минимайзеры: (" << std::setprecision(4) << x << ", " << y << ", " << z << ")";
        fout << " (функция была вычислена " << counter << " раз(а), её градиент - " << gradientCounter / 2 << " раз(а))"<< std::endl;
        fout << "• Корень полинома: (" << x << " + " << y << "i)" << std::endl;
        fout << std::setprecision(15);
    }
};
int function::counter = 0, function::gradientCounter = 0;


std::vector<double> CoordinateDescentBody(double x, double y, double epsilon, double step_x, double step_y,
                                          std::optional<double> x1 = std::nullopt, std::optional<double> y1 = std::nullopt, 
                                          std::optional<double> x2 = std::nullopt, std::optional<double> y2 = std::nullopt) {
    //Основная реализация метода покоординатного спуска
    int i = 0;
    fout << "┌────────┬──────────────────────┬──────────────────────┬──────────────────────┬───────┬──────────────────────┬──────────────────────┬─────────────────────┐\n"
            "│   i    │         x(i)         │         y(i)         │        F(x,y)        │ Коорд │        Шаг(х)        │        Шаг(у)        │    Трудоемкость     │\n"
            "├────────┼──────────────────────┼──────────────────────┼──────────────────────┼───────┼──────────────────────┼──────────────────────┼─────────────────────┤\n"
            ;
    double func = function::minimizer(x, y, x1, y1, x2, y2);
    double denom = 2, old_step_x = step_x, old_step_y= step_y;
    std::string type;
    fout << std::left;
    fout << "│ " << std::setw(7) << i << "│ " << std::setw(21)<< x << "│ " << std::setw(21) << y;
    fout << "│ " << std::setw(21)<< func;
    while (abs(func) >= epsilon) { //Критерий остановки
        bool flag1 = false, flag2 = false;
        type = "x";
        if (!function::stepAdvance(x, y, func, step_x, 0, step_x, step_y, i, type, x1, y1, x2, y2) &&
            !function::stepAdvance(x, y, func, -step_x, 0, -step_x, step_y, i, type, x1, y1, x2, y2)){
            step_x /= denom;
            flag1 = !flag1;
        }
        step_x /= 2;
        type = "y";
        if (!function::stepAdvance(x, y, func, 0, step_y, step_x, step_y, i, type, x1, y1, x2, y2) &&
            !function::stepAdvance(x, y, func, 0, -step_y, step_x, -step_y, i, type, x1, y1, x2, y2)) {
            step_y /= denom;
            flag2 = !flag2;
        }
        if (flag1 && flag2) {
            denom *= 2;
        }
        step_y /= 2;
        old_step_x = step_x; old_step_y = step_y;
    }
    fout << "│ " << std::setw(6) << type;
    fout << "│ " << std::setw(21)<< step_x << "│ " << std::setw(21) << step_y;
    fout << "│ " << std::setw(20)<< (i + 1) << "│ " << std::endl;
    fout << "└────────┴──────────────────────┴──────────────────────┴──────────────────────┴───────┴──────────────────────┴──────────────────────┴─────────────────────┘" << std::endl;
    std::vector<double> res;
    res.push_back(x); res.push_back(y);
    function::printResult(x, y, func);
    function::counter = 0; //Обнуляем счетчик
    fout << std::endl;
    return res;
}

std::vector<double> CrushingStepGradientBody(double x, double y, double epsilon, double alpha, double lambda,
                                             std::optional<double> x1 = std::nullopt, std::optional<double> y1 = std::nullopt,
                                             std::optional<double> x2 = std::nullopt, std::optional<double> y2 = std::nullopt) {
    //Основная реализация градиентного метода с дроблением шага
    int i = 0;
    fout << "┌────────┬──────────────────────┬──────────────────────┬──────────────────────┬──────────────────────┬──────────────────────┐\n"
            "│   i    │         x(i)         │         y(i)         │        F(x,y)        │       F′х(x,y)       │       F′у(x,y)       │\n"
            "├────────┼──────────────────────┼──────────────────────┼──────────────────────┼──────────────────────┼──────────────────────┤\n"
            ;
    double Func = function::minimizer(x, y, x1, y1, x2, y2);
    std::pair<double, double> Grad = std::make_pair(function::xDerivative(x, y, x1, y1, x2, y2),
                                                    function::yDerivative(x, y, x1, y1, x2, y2));
    double new_Func = function::minimizer(x - alpha * Grad.first, y - alpha * Grad.second, x1, y1, x2, y2);
    double squareNorm = pow(Grad.first, 2) + pow(Grad.second, 2);
    double z; // здесь будут значения ЦФ
    fout << std::left;
    fout << "│ " << std::setw(7) << i << "│ " << std::setw(21)<< x << "│ " << std::setw(21) << y;
    fout << "│ " << std::setw(21)<< Func;
    fout << "│ " << std::setw(21)<< Grad.first << "│ " << std::setw(21)<< Grad.second << "│ " << std::endl;
    i++;
    while (sqrt(squareNorm) >= epsilon) { //Критерий остановки
        while (new_Func - Func > - alpha * epsilon * squareNorm) {
            alpha *= lambda;
//            Func = function::minimizer(x, y, x1, y1, x2, y2);
            new_Func = function::minimizer(x - alpha * Grad.first, y - alpha * Grad.second, x1, y1, x2, y2);
        }

        x -= alpha * Grad.first;
        y -= alpha * Grad.second;
        Grad.first = function::xDerivative(x, y, x1, y1, x2, y2), Grad.second = function::yDerivative(x, y, x1, y1, x2, y2);
        z = function::minimizer(x, y, x1, y1, x2, y2); // ЦФ
        fout << "│ " << std::setw(7) << i << "│ " << std::setw(21)<< x << "│ " << std::setw(21) << y;
        fout << "│ " << std::setw(21)<< z;
        fout << "│ " << std::setw(21)<< Grad.first << "│ " << std::setw(21)<< Grad.second << "│ " << std::endl;
        i++;
        squareNorm = pow(Grad.first, 2) + pow(Grad.second, 2);
    }
    fout << "└────────┴──────────────────────┴──────────────────────┴──────────────────────┴──────────────────────┴──────────────────────┘" << std::endl;
    std::vector<double> res;
    res.push_back(x); res.push_back(y);
    function::printResult(x, y, z);
    function::counter = 0, function::gradientCounter = 0; //Обнуляем счетчики
    fout << std::endl;
    return res;
}

std::vector<double> ConstantStepGradientBody(double x, double y, double epsilon, double alpha,
                                             std::optional<double> x1 = std::nullopt, std::optional<double> y1 = std::nullopt,
                                             std::optional<double> x2 = std::nullopt, std::optional<double> y2 = std::nullopt) {
    //Основная реализация градиентного метода с постоянным шагом
    int i = 0;
    fout << "┌────────┬──────────────────────┬──────────────────────┬──────────────────────┬──────────────────────┬──────────────────────┐\n"
            "│   i    │         x(i)         │         y(i)         │        F(x,y)        │       F′х(x,y)       │       F′у(x,y)       │\n"
            "├────────┼──────────────────────┼──────────────────────┼──────────────────────┼──────────────────────┼──────────────────────┤\n"
            ;
    double Func = function::minimizer(x, y, x1, y1, x2, y2); // ЦФ
    std::pair<double, double> Grad = std::make_pair(function::xDerivative(x, y, x1, y1, x2, y2),
                                                    function::yDerivative(x, y, x1, y1, x2, y2));
    double norm = sqrt(pow(Grad.first, 2) + pow(Grad.second, 2));
    fout << std::left;
    fout << "│ " << std::setw(7) << i << "│ " << std::setw(21)<< x << "│ " << std::setw(21) << y;
    fout << "│ " << std::setw(21)<< Func;
    fout << "│ " << std::setw(21)<< Grad.first << "│ " << std::setw(21)<< Grad.second << "│ " << std::endl;
    i++;
    while (norm >= epsilon) { //Критерий остановки
        x -= alpha * Grad.first;
        y -= alpha * Grad.second;
        Func = function::minimizer(x, y, x1, y1, x2, y2); // ЦФ
        Grad.first = function::xDerivative(x, y, x1, y1, x2, y2), Grad.second = function::yDerivative(x, y, x1, y1, x2, y2);
        fout << "│ " << std::setw(7) << i << "│ " << std::setw(21)<< x << "│ " << std::setw(21) << y;
        fout << "│ " << std::setw(21)<< Func;
        fout << "│ " << std::setw(21)<< Grad.first << "│ " << std::setw(21)<< Grad.second << "│ " << std::endl;
        i++;
        norm = sqrt(pow(Grad.first, 2) + pow(Grad.second, 2));
    }
    fout << "└────────┴──────────────────────┴──────────────────────┴──────────────────────┴──────────────────────┴──────────────────────┘" << std::endl;
    std::vector<double> res;
    res.push_back(x); res.push_back(y);
    function::printResult(x, y, Func);
    function::counter = 0, function::gradientCounter = 0; //Обнуляем счетчики
    fout << std::endl;
    return res;
}

std::vector<double> PredefinedStepGradientBody(double x, double y, double epsilon, double k,
                                             std::optional<double> x1 = std::nullopt, std::optional<double> y1 = std::nullopt,
                                             std::optional<double> x2 = std::nullopt, std::optional<double> y2 = std::nullopt) {
    //Основная реализация градиентного спуска с заранее заданным шагом
    int i = 0;
    fout << "┌────────┬──────────────────────┬──────────────────────┬──────────────────────┬──────────────────────┬──────────────────────┐\n"
            "│   i    │         x(i)         │         y(i)         │        F(x,y)        │       F′х(x,y)       │       F′у(x,y)       │\n"
            "├────────┼──────────────────────┼──────────────────────┼──────────────────────┼──────────────────────┼──────────────────────┤\n"
            ;
    double Func = function::minimizer(x, y, x1, y1, x2, y2);
    std::pair<double, double> Grad = std::make_pair(function::xDerivative(x, y, x1, y1, x2, y2),
                                                    function::yDerivative(x, y, x1, y1, x2, y2));
    double norm = sqrt(pow(Grad.first, 2) + pow(Grad.second, 2));
    fout << std::left;
    fout << "│ " << std::setw(7) << i << "│ " << std::setw(21)<< x << "│ " << std::setw(21) << y;
    fout << "│ " << std::setw(21)<< Func;
    fout << "│ " << std::setw(21)<< Grad.first << "│ " << std::setw(21)<< Grad.second << "│ " << std::endl;
    i++;
    while (norm >= epsilon) { //Критерий остановки
        x -= Grad.first / k;
        y -= Grad.second / k;
        Func = function::minimizer(x, y, x1, y1, x2, y2);
        Grad.first = function::xDerivative(x, y, x1, y1, x2, y2), Grad.second = function::yDerivative(x, y, x1, y1, x2, y2);
        fout << "│ " << std::setw(7) << i << "│ " << std::setw(21)<< x << "│ " << std::setw(21) << y;
        fout << "│ " << std::setw(21)<< Func;
        fout << "│ " << std::setw(21)<< Grad.first << "│ " << std::setw(21)<< Grad.second << "│ " << std::endl;
        i++; k++;
        norm = sqrt(pow(Grad.first, 2) + pow(Grad.second, 2));
    }
    fout << "└────────┴──────────────────────┴──────────────────────┴──────────────────────┴──────────────────────┴──────────────────────┘" << std::endl;
    std::vector<double> res;
    res.push_back(x); res.push_back(y);
    function::printResult(x, y, Func);
    function::counter = 0, function::gradientCounter = 0; //Обнуляем счетчики
    fout << std::endl;
    return res;
}

std::vector<double> FastestGradientBody(double x, double y, double epsilon,
                                        double left, double right,
                                        double alphaLeft, double alphaRight,
                                        std::optional<double> x1 = std::nullopt, std::optional<double> y1 = std::nullopt,
                                        std::optional<double> x2 = std::nullopt, std::optional<double> y2 = std::nullopt) {
    //Основная реализация метода наискорейшего градиентного спуска
    int i = 0;
    fout << "┌────────┬──────────────────────┬──────────────────────┬──────────────────────┬──────────────────────┬──────────────────────┐\n"
            "│   i    │         x(i)         │         y(i)         │        F(x,y)        │       F′х(x,y)       │       F′у(x,y)       │\n"
            "├────────┼──────────────────────┼──────────────────────┼──────────────────────┼──────────────────────┼──────────────────────┤\n"
            ;
    double Func = function::minimizer(x, y, x1, y1, x2, y2);
    std::pair<double, double> Grad = std::make_pair(function::xDerivative(x, y, x1, y1, x2, y2),
                                                    function::yDerivative(x, y, x1, y1, x2, y2));
    double alpha;
    double norm = sqrt(pow(Grad.first, 2) + pow(Grad.second, 2));
    fout << std::left;
    fout << "│ " << std::setw(7) << i << "│ " << std::setw(21)<< x << "│ " << std::setw(21) << y;
    fout << "│ " << std::setw(21)<< Func;
    fout << "│ " << std::setw(21)<< Grad.first << "│ " << std::setw(21)<< Grad.second << "│ " << std::endl;
    i++;
    while (norm >= epsilon) { //Критерий остановки
        
        int k = ceil((right - left) / epsilon);
        double alphaMin = alphaLeft, FuncAlphaMin;
        for (int j = 0; j < k; j++) { //Критерий остановки
            double ai = alphaLeft + (alphaRight - alphaLeft) * (i + 1) / k;
            double FuncAlphai = function::alphaFunction(x, y, ai, x1, y1, x2, y2);
            if (FuncAlphai < FuncAlphaMin) {
                alphaMin = ai; FuncAlphaMin = FuncAlphai;
            } else
                break;
        }

        alpha = alphaMin;

        x -= alpha * Grad.first;
        y -= alpha * Grad.second;
        Func = function::minimizer(x, y, x1, y1, x2, y2);
        Grad.first = function::xDerivative(x, y, x1, y1, x2, y2), Grad.second = function::yDerivative(x, y, x1, y1, x2, y2);
        fout << "│ " << std::setw(7) << i << "│ " << std::setw(21)<< x << "│ " << std::setw(21) << y;
        fout << "│ " << std::setw(21)<< Func;
        fout << "│ " << std::setw(21)<< Grad.first << "│ " << std::setw(21)<< Grad.second << "│ " << std::endl;
        i++;
        norm = sqrt(pow(Grad.first, 2) + pow(Grad.second, 2));
    }
    fout << "└────────┴──────────────────────┴──────────────────────┴──────────────────────┴──────────────────────┴──────────────────────┘" << std::endl;
    std::vector<double> res;
    res.push_back(x); res.push_back(y);
    function::printResult(x, y, Func);
    function::counter = 0, function::gradientCounter = 0; //Обнуляем счетчики
    fout << std::endl;
    return res;
}

void CoordinateDescent(double x, double y, double epsilon, double step_x, double step_y) {
    //Метод покоординатного спуска
    fout << "МЕТОД ПОКООРДИНАТНОГО СПУСКА" << std::endl;
    fout << std::setprecision(15);
    //Находим минимум Cmin, делим на (z - Cmin) * (_z_ - _Сmin_) наш полином и дальше находим минимум нового полинома
    std::vector<double> root1 = CoordinateDescentBody(x, y, epsilon, step_x, step_y);
    std::vector<double> root2 = CoordinateDescentBody(x, y, epsilon, step_x, step_y, root1[0], root1[1]);
    CoordinateDescentBody(x, y, epsilon, step_x, step_y, root1[0], root1[1], root2[0], root2[1]);
}

void CrushingStepGradient(double x, double y, double epsilon, double alpha, double lambda) {
    //Градиентный метод с дроблением шага
    fout << "ГРАДИЕНТНЫЙ МЕТОД С ДРОБЛЕНИЕМ ШАГА" << std::endl;
    fout << std::setprecision(15);
    //Находим минимум Cmin, делим на (z - Cmin) * (_z_ - _Сmin_) наш полином и дальше находим минимум нового полинома
    std::vector<double> root1 = CrushingStepGradientBody(x, y, epsilon, alpha, lambda);
    std::vector<double> root2 = CrushingStepGradientBody(x, y, epsilon, alpha, lambda, root1[0], root1[1]);
    CrushingStepGradientBody(x, y, epsilon, alpha, lambda, root1[0], root1[1], root2[0], root2[1]);
}

void ConstantStepGradient(double x, double y, double epsilon, double alpha[3]) {
    //Градиентный метод с постоянным шагом
    fout << "ГРАДИЕНТНЫЙ МЕТОД С ПОСТОЯННЫМ ШАГОМ" << std::endl;
    fout << std::setprecision(15);
    //Находим минимум Cmin, делим на (z - Cmin) * (_z_ - _Сmin_) наш полином и дальше находим минимум нового полинома
    std::vector<double> root1 = ConstantStepGradientBody(x, y, epsilon, alpha[0]);
    std::vector<double> root2 = ConstantStepGradientBody(x, y, epsilon, alpha[1], root1[0], root1[1]);
    ConstantStepGradientBody(x, y, epsilon, alpha[2], root1[0], root1[1], root2[0], root2[1]);
}

void PredefinedStepGradient(double x, double y, double epsilon, double k[3]) {
    //Градиентный спуск с заранее заданным шагом
    fout << "ГРАДИЕНТНЫЙ СПУСК С ЗАРАНЕЕ ЗАДАННЫМ ШАГОМ" << std::endl;
    fout << std::setprecision(15);
    //Находим минимум Cmin, делим на (z - Cmin) * (_z_ - _Сmin_) наш полином и дальше находим минимум нового полинома
    std::vector<double> root1 = PredefinedStepGradientBody(x, y, epsilon, k[0]);
    std::vector<double> root2 = PredefinedStepGradientBody(x, y, epsilon, k[1], root1[0], root1[1]);
    PredefinedStepGradientBody(x, y, epsilon, k[2], root1[0], root1[1], root2[0], root2[1]);
}

void FastestGradient(double x, double y, double epsilon, double left, double right) {
    //МНГС
    fout << "МЕТОД НАИСКОРЕЙШЕГО ГРАДИЕНТНОГО СПУСКА" << std::endl;
    fout << std::setprecision(15);
    //Находим минимум Cmin, делим на (z - Cmin) * (_z_ - _Сmin_) наш полином и дальше находим минимум нового полинома
    std::vector<double> root1 = FastestGradientBody(x, y, epsilon, left, right, 0.0001, 0.2);
    std::vector<double> root2 = FastestGradientBody(0.1, 6, epsilon, left, right, 0.01, 0.2,root1[0], root1[1]);
    FastestGradientBody(8, 8, epsilon, left, right, 0.1, 0.2,root1[0], root1[1], root2[0], root2[1]);
}

int main() {

    CoordinateDescent(4, 4, 0.00001, 0.8, 0.9);
    CrushingStepGradient(4, 0.1, 0.00001, 1, 0.5);
    double alphaArray[3] = {0.0001, 0.01, 0.1};
    ConstantStepGradient(4, 1, 0.00001, alphaArray);
    double kArray[3] = {10000, 100, 1};
    PredefinedStepGradient(0, 0, 0.00001, kArray);
//    FastestGradient(6, 0.1, 0.00001, -1, 8);
    fout.close();
    return 0;
}
