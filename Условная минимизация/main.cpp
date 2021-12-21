#include <iostream>
#include <fstream>
//#include <complex>
#include <cmath>

#define PI 3.14159265

//Методы оптимизации: условная минимизация
//Функция (по z): (z - 7) * (z - 7i) * (z - 7 - 6i) = z^3 - (14 + 13i) * z^2 + (7 + 140i) * z + (294 - 343i)
//Функция (по x, y): x^3 + 3i*x^2*y - (14 + 13i)x^2 - 3xy^2 + (26 - 28i)xy + (7 + 140i)x - iy^3 + (14 + 13i)y^2 - (140 - 7i)y + (294 - 343i)
//Минимайзер (по х, у): (x^2 - 14x + y^2 + 49) * (x^2 + y^2 - 14y + 49) * (x^2 - 14x + y^2 - 12y + 85);
//Условие (по х, у): x^2 + y^2 - R^2 <= 0;

std::ofstream fout("output.txt");
const double R = 4; //Радиус окружности, внутри которой ищем точку минимума

struct function {
    static int functionCounter, constrainCounter, gradientCounter, constrainGradientCounter;

    static double minimizer(double x, double y) {
        //Возвращает значение функции |P(x, y)|^2 в точке (x, y)
        functionCounter++; //Счетчик вычислений функций
        return (pow(x, 2) - 14 * x + pow(y, 2) + 49) * (pow(x, 2) + pow(y, 2) - 14 * y + 49) * (pow(x, 2) - 14 * x + pow(y, 2) - 12 * y + 85);
    }

    static double xDerivative(double x, double y) {
        //Возвращает значение производной функции |P(x, y)|^2 по х в точке (x, y)
        gradientCounter++; //Счетчик вычислений производной
        return 6*pow(x, 5) - 140*pow(x, 4) + 12*pow(x, 3)*pow(y, 2) - 104*pow(x, 3)*y + 1516*pow(x, 3) - 168*pow(x, 2)*pow(y, 2) + 1680*pow(x, 2)*y - 9744*pow(x, 2) + 6*x*pow(y, 4) - 104*x*pow(y, 3) + 1460*x*pow(y, 2) - 11592*x*y + 40670*x - 28*pow(y, 4) + 560*pow(y, 3) - 5600*pow(y, 2) + 34496*y - 91924;
    }

    static double yDerivative(double x, double y) {
        //Возвращает значение производной функции |P(x, y)|^2 по y в точке (x, y)
        gradientCounter++; //Счетчик вычислений функций
        return 6*pow(x, 4)*y - 26*pow(x, 4) - 112*pow(x, 3)*y + 560*pow(x, 3) + 12*pow(x, 2)*pow(y, 3) - 156*pow(x, 2)*pow(y, 2) + 1460*pow(x, 2)*y - 5796*pow(x, 2) - 112*x*pow(y, 3) + 1680*x*pow(y, 2) - 11200*x*y + 34496*x + 6*pow(y, 5) - 130*pow(y, 4) + 1404*pow(y, 3) - 9156*pow(y, 2) + 37926*y - 87122;
    }

    static double constrain(double x, double y) {
        //Возвращает значение условия g(x, y) в точке (x, y)
        constrainCounter++; //Счетчик вычислений функций
        return pow(x, 2) + pow(y, 2) - pow(R, 2);
    }

    static double xConstrainDerivative(double x, double y) {
        //Возвращает значение производной функции |P(x, y)|^2 по х в точке (x, y)
        constrainGradientCounter++; //Счетчик вычислений производной
        return 2 * x;
    }

    static double yConstrainDerivative(double x, double y) {
        //Возвращает значение производной функции |P(x, y)|^2 по y в точке (x, y)
        constrainGradientCounter++; //Счетчик вычислений функций
        return 2 * y;
    }

    static double norm(std::pair<double, double> vector) {
        //Возвращает значение нормы вектора {x, y}
        return sqrt(pow(vector.first, 2) + pow(vector.second, 2));
    }

    static double dotProduct(std::pair<double, double> vector1, std::pair<double, double> vector2) {
        //Возвращает скалярное произведение двух векторов {x1, y1} и {x2, y2}
        return vector1.first * vector2.first + vector1.second * vector2.second;
    }

    static void printResult(double x, double y, double z){
        fout << "• Минимум в точке: (" << x << ", " << y << "), и он равен: " << z << std::endl;
        fout << "• Приближенные минимайзеры: (" << std::setprecision(4) << x << ", " << y << ", " << z << ")";
        fout << " (функция f была вычислена " << functionCounter << " раз(а), её градиент - " << gradientCounter / 2;
        fout << " раз(а); функция g была вычислена " << constrainCounter << " раз(а), её градиент - " << constrainGradientCounter / 2 << std::endl;
        fout << "• g(x, y): x² + y² - " << pow(R, 2) <<" ≤ 0" << std::endl;
        fout << "• Угол в критерии остановки: 4.24728938°" << std::endl;
        fout << std::setprecision(6);
    }
};
int function::functionCounter = 0, function::constrainCounter = 0,
    function::gradientCounter = 0, function::constrainGradientCounter = 0;

void ConstrainedMinimization(double x, double y, double epsilon, double alpha) {
    //Условная минимизация
    fout << "УСЛОВНАЯ МИНИМИЗАЦИЯ" << std::endl;
    fout << std::setprecision(6);
    int i = 0;
    fout << "┌────────┬──────────────┬──────────────┬──────────────┬───────────────────────────────┬──────────────┬───────────────────────────────┬──────────────┬───────────────────────────────┬──────────────┐\n"
            "│   i    │     x(i)     │     y(i)     │    f(x,y)    │            ∇f(x,y)            │    g(x,y)    │            ∇g(x,y)            │  ∠(-∇f, ∇g)  │            pr(-∇f)            │ Множитель  α │\n"
            "├────────┼──────────────┼──────────────┼──────────────┼───────────────────────────────┼──────────────┼───────────────────────────────┼──────────────┼───────────────────────────────┼──────────────┤\n"
            ;
    double Func = function::minimizer(x, y);
    std::pair<double, double> Grad = std::make_pair(function::xDerivative(x, y),
                                                    function::yDerivative(x, y));
    std::pair<double, double> antiGrad = std::make_pair(-Grad.first, -Grad.second);
    double Constrain = function::constrain(x, y);
    std::pair<double, double> constrainGrad = std::make_pair(function::xConstrainDerivative(x, y),
                                                             function::yConstrainDerivative(x, y));
    double currentCos = function::dotProduct(antiGrad, constrainGrad) / (function::norm(Grad) * function::norm(constrainGrad));
    fout << std::left;
    fout << "│ " << std::setw(7) << i << "│ " << std::setw(13) << x << "│ " << std::setw(13)  << y;
    fout << "│ " << std::setw(13) << Func;
    fout << "│ x: " << std::setw(11)  << Grad.first << "y: " << std::setw(13) << Grad.second;
    fout << "│ " << std::setw(13) << Constrain;
    fout << "│ x: " << std::setw(11)  << constrainGrad.first << "y: " << std::setw(13) << constrainGrad.second;
    fout << "│ " << std::setw(13)  << acos(currentCos) * 180 / PI;
    fout << "│ x: " << std::setw(11)  << "-" << "y: " << std::setw(13) << "-";
    fout << "│ " << std::setw(13) << alpha << "│ "<< std::endl;
    i++;

    while (abs(Constrain) > epsilon) {
        //Спуск до границы
        Constrain = function::constrain(x + alpha * antiGrad.first, y + alpha * antiGrad.second);
        if (Constrain > 0) {
            alpha /= 10;
            continue;
        }
        x += alpha * antiGrad.first;
        y += alpha * antiGrad.second;
        Func = function::minimizer(x, y);
        Grad.first = function::xDerivative(x, y), Grad.second = function::yDerivative(x, y);
        antiGrad.first = -Grad.first, antiGrad.second = -Grad.second;
        constrainGrad.first = function::xConstrainDerivative(x, y);
        constrainGrad.second = function::yConstrainDerivative(x, y);
        currentCos = function::dotProduct(antiGrad, constrainGrad) / (function::norm(Grad) * function::norm(constrainGrad));
        fout << "│ " << std::setw(7) << i << "│ " << std::setw(13) << x << "│ " << std::setw(13)  << y;
        fout << "│ " << std::setw(13) << Func;
        fout << "│ x: " << std::setw(11)  << Grad.first << "y: " << std::setw(13) << Grad.second;
        fout << "│ " << std::setw(13) << Constrain;
        fout << "│ x: " << std::setw(11)  << constrainGrad.first << "y: " << std::setw(13)<< constrainGrad.second;
        fout << "│ " << std::setw(13)  << acos(currentCos) * 180 / PI;
        fout << "│ x: " << std::setw(11)  << "-" << "y: " << std::setw(13) << "-";
        fout << "│ " << std::setw(13) << alpha << "│ "<< std::endl;
        i++;
    }

//    fout << "├────────┼──────────────────────┼──────────────────────┼──────────────────────┼───────────────────────────────────────────────┼──────────────────────┼───────────────────────────────────────────────┼──────────────────────┼───────────────────────────────────────────────┤\n";
    const double cos = 0.99725369; //чуть меньше, чем cos5°
    std::pair<double, double> projection = std::make_pair(constrainGrad.first * (function::dotProduct(antiGrad, constrainGrad) / function::dotProduct(constrainGrad, constrainGrad)),
                                                          constrainGrad.second * (function::dotProduct(antiGrad, constrainGrad) / function::dotProduct(constrainGrad, constrainGrad)));
    while (currentCos < cos) {
        //Скольжение по границе
        //std::cout << std::setprecision(6) << (function::dotProduct(antiGrad, constrainGrad) / (function::norm(Grad) * function::norm(constrainGrad))) << std::endl;
        x += alpha * (antiGrad.first - projection.first);
        y += alpha * (antiGrad.second - projection.second);

        Func = function::minimizer(x, y);
        Grad.first = function::xDerivative(x, y), Grad.second = function::yDerivative(x, y);
        antiGrad.first = -Grad.first, antiGrad.second = -Grad.second;
        Constrain = function::constrain(x, y);
        constrainGrad.first = function::xConstrainDerivative(x, y);
        constrainGrad.second = function::yConstrainDerivative(x, y);
        currentCos = function::dotProduct(antiGrad, constrainGrad) / (function::norm(Grad) * function::norm(constrainGrad));
        projection.first = constrainGrad.first * (function::dotProduct(antiGrad, constrainGrad) / function::dotProduct(constrainGrad, constrainGrad));
        projection.second = constrainGrad.second * (function::dotProduct(antiGrad, constrainGrad) / function::dotProduct(constrainGrad, constrainGrad));
        fout << "│ " << std::setw(7) << i << "│ " << std::setw(13) << x << "│ " << std::setw(13)  << y;
        fout << "│ " << std::setw(13) << Func;
        fout << "│ x: " << std::setw(11)  << Grad.first << "y: " << std::setw(13) << Grad.second;
        fout << "│ " << std::setw(13) << Constrain;
        fout << "│ x: " << std::setw(11)  << constrainGrad.first << "y: " << std::setw(13) << constrainGrad.second;
        fout << "│ " << std::setw(13)  << acos(currentCos) * 180 / PI;
        fout << "│ x: " << std::setw(11)  << projection.first << "y: " << std::setw(13) << projection.second;
        fout << "│ " << std::setw(13) << alpha << "│ "<< std::endl;
        i++;
    }
    fout << "└────────┴──────────────┴──────────────┴──────────────┴───────────────────────────────┴──────────────┴───────────────────────────────┴──────────────┴───────────────────────────────┴──────────────┘" << std::endl;
    function::printResult(x, y, Func);
    function::functionCounter = 0, function::constrainCounter = 0,
    function::gradientCounter = 0, function::constrainGradientCounter = 0; //Обнуляем счетчики
}


int main() {
    ConstrainedMinimization(0.1, 0.1, 0.00001, 0.0001);
    return 0;
}
