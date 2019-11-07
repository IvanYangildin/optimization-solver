#ifndef PARAMSDIALOG_H
#define PARAMSDIALOG_H

#include <QWidget>
#include "IProblem.h"
#include "ISolver.h"
#include <QDialog>
#include <QString>
#include "IVector.h"

namespace Ui {
    class ParamsDialog;
}

class ParamsDialog : public QDialog
{
    Q_OBJECT

public:
    explicit ParamsDialog(QWidget *parent, IProblem* problem, ISolver* solver);
    ~ParamsDialog();

private slots:
    void on_solve_button_released();

    void on_by_args_button_clicked();

    void on_by_params_button_clicked();

    void on_solve_button_clicked();

private:
    QString mergeString(size_t dim_args, size_t dim_params, const QString * epsilon,
                         int flag,  const QString * start_point, const QString * vector, const QString * begin, const QString * end);

    Ui::ParamsDialog *ui;
    IProblem* m_problem;
    ISolver* m_solver;
};

#endif // PARAMSDIALOG_H
