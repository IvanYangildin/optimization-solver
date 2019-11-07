#include "paramsdialog.h"
#include "ui_paramsdialog.h"
#include <QMessageBox>
#include "ILog.h"
#include "ICompact.h"
#include "qplot.h"
#include "ISet.h"
#include <math.h>
#include "error.h"

ParamsDialog::ParamsDialog(QWidget *parent, IProblem* problem, ISolver* solver) :
    QDialog(parent),
    ui(new Ui::ParamsDialog)
{
    m_problem = problem;
    m_solver = solver;
    ui->setupUi(this);

    size_t dim_left = -1, dim_right = -1;
    if (m_problem->getArgsDim(dim_left) != ErrorEnum::ERR_OK) {
        QMessageBox::warning(this, tr("error"), tr("can't get dim_right"));
    }
    if (m_problem->getParamsDim(dim_right) != ErrorEnum::ERR_OK) {
        QMessageBox::warning(this, tr("error"), tr("can't get dim_right"));
    }
    ui->dim_left_title->setText("Dimension is: " + QString("%1").arg(dim_left));
    ui->dim_right_title->setText("Dimension is: " + QString("%1").arg(dim_right));
}

ParamsDialog::~ParamsDialog()
{
    delete ui;
}

QString ParamsDialog::mergeString(size_t dim_args, size_t dim_params, const QString * epsilon,
                                  int flag,  const QString * start_point, const QString * vector, const QString * begin, const QString * end) {
    // TODO
    QString comma = ", ";
    QString res;

    res = QString("%1").arg(dim_args);
    res += comma;
    res += QString("%1").arg(dim_params);
    res += comma;
    res += epsilon;
    res += comma;
    res += QString("%1").arg(flag);
    res += comma;


    if (flag) {
        res += start_point;
        res += comma;
        res += vector;
    }
    else {
        res += vector;
        res += comma;
        res += start_point;
    }

    res += comma;
    res += begin;
    res += comma;
    res += end;

    return res;
}

IVector* StrToVec(const QString *str) {
    QStringList coordsList = str->split(',');
    double *vecD = new double[coordsList.size()];

    for (int i = 0; i < coordsList.size(); i++) {
        bool res_convert = false;
        double elem = coordsList.at(i).toDouble(&res_convert);
        if (res_convert == false) {
            ILog::report("StringToVector: can't convert str to double");
            return NULL;
        }
        vecD[i] = elem;
    }

    IVector *vecParams = IVector::createVector(coordsList.size(), vecD);
    return vecParams;
}

void ParamsDialog::on_solve_button_released()
{
    if (!m_problem) {
        QMessageBox::warning(this, tr("error"), tr("problem is not exist!"));
        ILog::report("problem is not exist!");
        return;
    }

    if (!m_solver) {
        QMessageBox::warning(this, tr("error"), tr("solver is not exist!"));
        ILog::report("solver is not exist!");
        return;
    }


    bool isByArgs;
    if (ui->by_args_button->isChecked()) {
        isByArgs = true;
    }
    else if (ui->by_params_button->isChecked()) {
        isByArgs = false;
    }
    else {
        QMessageBox::warning(this, tr("error"), tr("please choose radio button"));
        ILog::report("please choose radio button");
        return;
    }

    size_t dim_left, dim_right;

    if (m_problem->getArgsDim(dim_left) != ErrorEnum::ERR_OK) {
        QMessageBox::warning(this, tr("error"), tr("can't get dim_left"));
        ILog::report("can't get dim_left");
        return;
    }

    if (m_problem->getParamsDim(dim_right) != ErrorEnum::ERR_OK) {
        QMessageBox::warning(this, tr("error"), tr("can't get dim_right"));
        ILog::report("can't get dim_right");
        return;
    }

    size_t dim_args, dim_params;
    if (isByArgs) {
        dim_args = dim_left;
        dim_params = dim_right;
    }
    else {
        dim_args = dim_right;
        dim_params = dim_left;
    }

    const QString epsilon = ui->epsilon_editor->text();
    const QString start_point = ui->start_point_editor->text();
    const QString vector = ui->vector_editor->text();
    const QString begin = ui->begin_editor->text();
    const QString end = ui->end_editor->text();
    ISet* set = ISet::createSet(dim_left);
    IVector* beginV = StrToVec(&begin);
    IVector* endV = StrToVec(&end);
    IVector* start_pointV = StrToVec(&start_point);


    for (int i = 0; i < dim_left; i++) {
        double a, b, point;
        if(beginV->getCoord(i, a) != ERR_OK)
        {
            QMessageBox::warning(this, tr("error"), tr("can`t get coord of begin"));
            ILog::report("can`t get coord of begin");
            return;
        }
        if(endV->getCoord(i, b) != ERR_OK)
        {
            QMessageBox::warning(this, tr("error"), tr("can`t get coord of end"));
            ILog::report("can`t get coord of end");
            return;
        }
        if(start_pointV->getCoord(i, point) != ERR_OK)
        {
            QMessageBox::warning(this, tr("error"), tr("can`t get coord of end"));
            ILog::report("can`t get coord of end");
            return;
        }

        if (b < a) {
            QMessageBox::warning(this, tr("error"), tr("end[i] < begin[i]"));
            ILog::report("begin[i] > end[i]");
            return;
        }

        if (!(a <= point && point <= b)) {
            QMessageBox::warning(this, tr("error"), tr("start point is not contains in compact!"));
            ILog::report("start point is not contains in compact!");
            return;
        }

    }

    int flag = isByArgs ? 1 : 0;

    QString merged_string = mergeString(dim_args, dim_params, &epsilon, flag, &start_point, &vector, &begin, &end);

    if (m_solver->setProblem(m_problem) != ERR_OK) {
        QMessageBox::warning(this, tr("error"), tr("can`t set problem"));
        ILog::report("can't set problem");
        return;
    }
    else{
        //QMessageBox::information(this, tr("ok"), tr("set problem ok"));
    }

    if (m_solver->setParams(merged_string) != ERR_OK) {
        QMessageBox::warning(this, tr("error"), tr("can`t set params"));
        ILog::report("can't set params");
        return;
    }
    else{
        //QMessageBox::information(this, tr("ok"), tr("set params ok"));
    }

    if (m_solver->solve() != ERR_OK){
        QMessageBox::warning(this, tr("error"), tr("can't solve"));
        ILog::report("can't solve");
        return;
    }
    else{
        //QMessageBox::information(this, tr("ok"), tr("solve is done"));
    }

    double * solutionD = new (std::nothrow) double[dim_left];
    IVector * solution = IVector::createVector(dim_left, solutionD);
    m_solver->getSolution(solution);
    set->put(solution);

    IVector* v = IVector::createVector(dim_left, solutionD);
    ISet::IIterator* it = set->begin();

    do{
        if (set->getByIterator(it, v) != ERR_OK) {
            ILog::report("ERR: can't get by iterator in set");
            break;
        }
        it->next();
        QPlot* form = new QPlot(this);
        form->set(m_problem, v, StrToVec(&vector), beginV, endV, isByArgs);
        form->show();
    } while (!it->isEnd());

    delete set;
    delete beginV;
    delete endV;
    delete v;
    delete solutionD;
    delete solution;
}


void ParamsDialog::on_by_args_button_clicked()
{
    size_t dim_args, dim_params;
    if (m_problem->getArgsDim(dim_args) != ErrorEnum::ERR_OK) {
        QMessageBox::warning(this, tr("error"), tr("can't get dim_left"));
        ILog::report("can't get dim_left");
        return;
    }

    if (m_problem->getParamsDim(dim_params) != ErrorEnum::ERR_OK) {
        QMessageBox::warning(this, tr("error"), tr("can't get dim_right"));
        ILog::report("can't get dim_right");
        return;
    }

    QString text = "Dimension is: ";
    ui->dim_left_title->setText(text + QString("%1").arg(dim_args));
    ui->dim_right_title->setText(text + QString("%1").arg(dim_params));
    ui->title_left->setText(tr("Arguments"));
    ui->title_right->setText(tr("Parameters"));
}


void ParamsDialog::on_by_params_button_clicked()
{
    size_t dim_args, dim_params;
    if (m_problem->getArgsDim(dim_args) != ErrorEnum::ERR_OK) {
        QMessageBox::warning(this, tr("error"), tr("can't get dim_left"));
        ILog::report("can't get dim_left");
        return;
    }

    if (m_problem->getParamsDim(dim_params) != ErrorEnum::ERR_OK) {
        QMessageBox::warning(this, tr("error"), tr("can't get dim_right"));
        ILog::report("can't get dim_right");
        return;
    }

    QString text = "Dimension is: ";
    ui->dim_left_title->setText(text + QString("%1").arg(dim_params));
    ui->dim_right_title->setText(text + QString("%1").arg(dim_args));
    ui->title_left->setText(tr("Parameters"));
    ui->title_right->setText(tr("Arguments"));
}

void ParamsDialog::on_solve_button_clicked()
{

}
