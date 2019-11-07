#include "qplot.h"
#include "ui_qplot.h"
#include <QFileDialog>
#include <QFileInfo>
#include <QtSql>
#include <QLibrary>
#include <QMessageBox>
#include <float.h>
#include <limits.h>
#include "qcustomplot.h"
#include <iostream>

QPlot::QPlot(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::QPlot)
{
    ui->setupUi(this);
}

QPlot::~QPlot()
{
    delete _solution;
    delete _variables;
    delete _begin;
    delete _end;
    delete ui;
}

void QPlot::on_drawButton_clicked()
{
    _curr_axis = ui->spinAxis->text().toUInt();


    _a = ui->spinLeftBorder->text().toDouble();
    _b = ui->spinRightBorder->text().toDouble();

    double a_beg, b_end;

    if(_begin->getCoord(_curr_axis - 1, a_beg) != ERR_OK)
    {
        QMessageBox::critical(this, tr("Failed to draw plot"), tr("can`t get coord of begin"));
        return;
    }

    if(_end->getCoord(_curr_axis - 1, b_end) != ERR_OK)
    {
        QMessageBox::critical(this, tr("Failed to draw plot"), tr("can`t get coord of end"));
        return;
    }

    _a = std::max(_a, a_beg);
    _b = std::min(_b, b_end);
    double h = 0.01, ah = (_b - _a) / (INT_MAX - 2);
    size_t N = (_b - _a) / h + 2;

    if (_a > _b) {
        QMessageBox::critical(this, tr("Failed to draw plot"), tr("Bad borders."));
        return;
    }

    if (ah > h) {
        N = INT_MAX;
        h = ah;
    }

    QVector<double> x(N), y(N);
    QVector<double> sol_x(1), sol_y(1);
    IVector * iter = _solution->clone();
    if (!iter) {
        QMessageBox::critical(this, tr("Failed to draw plot"), tr("Failed to alloc memory."));
        return;
    }

    int last_error;


    double X = _a;
    for (int i = 0; i < N; i++)
    {
        iter->setCoord(_curr_axis - 1, X);

        x[i] = X;

        if (_isByArgs) {
            last_error = _problem->goalFunction(iter, _variables, y[i]);
        } else {
            last_error = _problem->goalFunction(_variables, iter, y[i]);
        }

        if (last_error != ERR_OK) {
            delete iter;
            QMessageBox::critical(this, tr("Failed to draw plot"), tr("Failed to get value of goal function."));
            return;
        }
        X += h;
    }

    delete iter;
    if (_solution->getCoord(_curr_axis - 1, sol_x[0]) != ERR_OK) {
        QMessageBox::critical(this, tr("Failed to draw plot"), tr("Failed to get value of solution point."));
        return;
    }

    if (_isByArgs) {
        last_error = _problem->goalFunction(_solution, _variables, sol_y[0]);
    } else {
        last_error = _problem->goalFunction(_variables, _solution, sol_y[0]);
    }
    if (last_error != ERR_OK) {
        QMessageBox::critical(this, tr("Failed to draw plot"), tr("Failed to get value of goal function in the solution point."));
        return;
    }

    ui->plot->clearGraphs();

    ui->plot->addGraph();

    //for(int i = 0; i < N; i++)
        //std::cout << i << ") x[" << i << "] = " << x[i] << "\ty[" << i << "] = " << y[i] << std::endl;

    ui->plot->graph(0)->setData(x, y);

    ui->plot->xAxis->setLabel(QString("x_%1").arg(_curr_axis));
    ui->plot->yAxis->setLabel("f");
    ui->plot->xAxis->setRange(_a, _b);

    double minY = y[0], maxY = y[0];
    for (int i = 1; i < N; ++i)
    {
        if (y[i] < minY) minY = y[i];
        else if (y[i] > maxY) maxY = y[i];
    }
    ui->plot->yAxis->setRange(minY, maxY);

    ui->plot->addGraph(ui->plot->xAxis, ui->plot->yAxis);
    ui->plot->graph(1)->setPen(QColor(50, 50, 50, 255));
    ui->plot->graph(1)->setLineStyle(QCPGraph::lsNone);
    ui->plot->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCross, 16));
    ui->plot->graph(1)->setData(sol_x, sol_y);

    ui->plot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);

    ui->plot->replot();
}


void QPlot::set(IProblem *problem, IVector *solution, IVector *variables,
                IVector *begin, IVector *end, bool isByArgs)
{
    _solution = solution->clone();
    _problem = problem;
    _variables = variables->clone();
    _begin = begin->clone();
    _end = end->clone();
    _isByArgs = isByArgs;
    unsigned int dim  = solution->getDim();
    ui->spinAxis->setMaximum(dim);
    ui->spinAxis->setMinimum(1);

    double min_margin, max_margin;
    _begin->getCoord(0, min_margin);
    _end->getCoord(0, max_margin);

    ui->spinLeftBorder->setMinimum(min_margin);
    ui->spinLeftBorder->setMaximum(max_margin);
    ui->spinRightBorder->setMinimum(min_margin);
    ui->spinRightBorder->setMaximum(max_margin);
    ui->spinLeftBorder->setValue(min_margin);
    ui->spinRightBorder->setValue(max_margin);
}

void QPlot::on_spinAxis_valueChanged(const QString &arg1)
{
    _curr_axis = ui->spinAxis->text().toUInt();
    double min_margin, max_margin;
    _begin->getCoord(_curr_axis, min_margin);
    _end->getCoord(_curr_axis, max_margin);

    ui->spinLeftBorder->setMinimum(min_margin);
    ui->spinLeftBorder->setMaximum(max_margin);
    ui->spinRightBorder->setMinimum(min_margin);
    ui->spinRightBorder->setMaximum(max_margin);

}
