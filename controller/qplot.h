#ifndef QPLOT_H
#define QPLOT_H

#include <QWidget>
#include <QDialog>
#include "IProblem.h"
#include "ISolver.h"
#include "ISet.h"
#include "ILog.h"

namespace Ui {
    class QPlot;
}

class QPlot : public QDialog
{
    Q_OBJECT

public:
    explicit QPlot(QWidget *parent);

    void set(IProblem *problem, IVector *solution, IVector *variables,
             IVector *begin, IVector *end, bool isByArgs);
    ~QPlot();

private slots:
    void on_drawButton_clicked();

    void on_spinAxis_valueChanged(const QString &arg1);

private:
    Ui::QPlot *ui;
    ISolver *_solver;
    IProblem *_problem;
    bool _isByArgs;
    IVector* _solution;
    IVector* _variables;
    IVector* _begin;
    IVector* _end;

    size_t _curr_axis;
    double _a, _b;
};

#endif // QPLOT_H
