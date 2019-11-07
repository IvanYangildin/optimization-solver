#include "controller.h"
#include "sqlconnectiondialog.h"
#include <QFileDialog>
#include <QFileInfo>
#include <QtSql>
#include <QMessageBox>
#include <float.h>
#include <limits.h>
#include "qplot.h"
#include "ILog.h"
#include "IBrocker.h"
#include "paramsdialog.h"

#include <QWidget>
#include <QCloseEvent>
Controller::Controller(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::Controller)
{
    ui->setupUi(this);
    windows = new QList<QWidget*>();
    m_problem = NULL;
    m_solver = NULL;
    m_problem_brocker = NULL;
    m_solver_brocker = NULL;
}

Controller::~Controller()
{
    if (m_problem_brocker)
        m_problem_brocker->release();
    if (m_solver_brocker)
        m_solver_brocker->release();
    delete ui;
}

QSqlError Controller::changeConnection(const QString &driver, const QString &dbName, const QString &host,
                             const QString &user, const QString &passwd, int port)
{
    static int cCount = 0;

    QSqlError err;
    QSqlDatabase db = QSqlDatabase::addDatabase(driver, QString("optimisation_solver%1").arg(++cCount));
    db.setDatabaseName(dbName);
    db.setHostName(host);
    db.setPort(port);
    if (!db.open(user, passwd)) {
        err = db.lastError();
        db = QSqlDatabase();
        QSqlDatabase::removeDatabase(QString("optimization_solver%1").arg(cCount));
        QMessageBox::warning(this, tr("Unable to open database"),tr("1"));
    } else if (select()){
        //QMessageBox::warning(this, tr("msg"),tr("DB has opened"));
        //insButton->setEnabled(true);
    } else {
        //insButton->setEnabled(false);
        QMessageBox::warning(this, tr("Unable to open database"),tr("3"));
        db.close();
        QSqlDatabase::removeDatabase(QString("optimization_solver%1").arg(cCount));
    }
    return err;
}

void Controller::changeConnection(){
    SqlConnectionDialog dialog(this);
    if (dialog.exec() != QDialog::Accepted)
        return;
    QSqlError err = changeConnection(dialog.driverName(), dialog.databaseName(), dialog.hostName(),
                                     dialog.userName(), dialog.password(), dialog.port());
    if (err.type() != QSqlError::NoError)
        QMessageBox::warning(this, tr("Unable to open database"), tr("An error occurred while "
                             "opening the connection: ") + err.text());

}

bool Controller::select()
{
    QStringList connectionNames = QSqlDatabase::connectionNames();// список всех подключений, возьмем последнее
    QSqlDatabase currentDB = QSqlDatabase::database(connectionNames.value(connectionNames.count() - 1));
    if (currentDB.tables().count() == 0) {
       QSqlQueryModel *m = new QSqlQueryModel(ui->table);
       m->setQuery(QSqlQuery("create table solvers(Name text, Description text, Path text);", currentDB));
       if (m->lastError().type() != QSqlError::NoError) {
           QMessageBox::warning(this, tr("Unable to open database"),tr("db3_1"));
           emit statusMessage(m->lastError().text());
           delete m;
           return false;
       }
       delete m;
    }
    QStringList tables = currentDB.tables();
    QSqlTableModel *model = new QSqlTableModel(ui->table, currentDB);
    model->setEditStrategy(QSqlTableModel::OnRowChange);
    model->setTable(currentDB.driver()->escapeIdentifier(tables.at(0), QSqlDriver::TableName));
    model->select();
    if (model->lastError().type() != QSqlError::NoError) {
        QMessageBox::warning(this, tr("Unable to open database"),tr("db3_2"));
        emit statusMessage(model->lastError().text());
        delete model;
        return false;
    }
    if (!checkTableModel(model->record())) {
        QMessageBox::warning(this, tr("Unable to open database"),tr("db3_3"));
        emit statusMessage(tr("Bad table."));
        delete model;
        return false;
    }
    ui->table->setModel(model);
    ui->table->setEditTriggers(QAbstractItemView::DoubleClicked|QAbstractItemView::EditKeyPressed);

    connect(ui->table->selectionModel(), SIGNAL(currentRowChanged(QModelIndex, QModelIndex)),
            this, SLOT(currentChanged()));

    if (tables.count() > 1)
        emit statusMessage(tr("Query OK. Count of tables is more than one. Using first table."));
    else
        emit statusMessage(tr("Query OK."));

    updateActions();
    return true;
}

void Controller::updateActions()
{
    bool enableIns = qobject_cast<QSqlTableModel *>(ui->table->model());
    bool enableDel = enableIns && ui->table->currentIndex().isValid();

    ui->action_Insert_Row->setEnabled(true);
    ui->action_Delete_Row->setEnabled(true);
}


bool Controller::checkTableModel(QSqlRecord record) {
    QString fieldNames[3] = {"Name", "Description", "Path"};
    if (record.count() != 3)
        return false;
    for (int i = 0; i < 3; i++){
        if (record.fieldName(i).compare(fieldNames[i]) != 0){
            return false;
        }
    }
    return true;
}

void Controller::insertRow()
{
    QSqlTableModel *model = qobject_cast<QSqlTableModel *>(ui->table->model());
    if (!model)
        return;

    QModelIndex insertIndex = ui->table->currentIndex();
    int row = insertIndex.row() == -1 ? 0 : insertIndex.row();
    model->insertRow(row);
    insertIndex = model->index(row, 0);
    ui->table->setCurrentIndex(insertIndex);
    ui->table->edit(insertIndex);
}

void Controller::deleteRow()
{
    QSqlTableModel *model = qobject_cast<QSqlTableModel *>(ui->table->model());
    if (!model)
        return;

    model->setEditStrategy(QSqlTableModel::OnManualSubmit);

    QModelIndexList currentSelection = ui->table->selectionModel()->selectedIndexes();
    for (int i = 0; i < currentSelection.count(); ++i) {
        if (currentSelection.at(i).column() != 0)
            continue;
        model->removeRow(currentSelection.at(i).row());
    }

    model->submitAll();
    model->setEditStrategy(QSqlTableModel::OnRowChange);

    updateActions();
}

void Controller::on_SolveButton_clicked()
{

}

bool Controller::loadBrocker(IBrocker::Type type, const QString &path) {
    switch (type)
    {
    case IBrocker::PROBLEM:
        if (m_problem_brocker) {
            m_problem_brocker->release();
            m_problem = NULL;
            m_problem_brocker = NULL;
        }
        break;
    case IBrocker::SOLVER:
      if (m_solver_brocker) {
          m_solver_brocker->release();
          m_solver = NULL;
          m_solver_brocker = NULL;
      }
      break;
    default:
        QMessageBox::warning(this, tr("error"), tr("wrong type argument"));
        ILog::report("wrong type argument");
        return false;
    }

    QLibrary lib(path);
    if (!lib.load()) {
        QMessageBox::warning(this, tr("error"), QString("can't load brocker") + path);
        ILog::report("can't load brocker");
        return false;
    }

    GetBrockerFunc getBrockerFunc = reinterpret_cast<GetBrockerFunc>(lib.resolve("getBrocker"));
    if (!getBrockerFunc) {
        QMessageBox::warning(this, tr("error"), tr("can't resolve brocker"));
        ILog::report("can't resolve brocker");
        return false;
    }

    IBrocker *brocker = reinterpret_cast<IBrocker *>(getBrockerFunc());
    void * interfaceImpl = NULL;
    if (brocker->canCastTo(type)) {
        interfaceImpl = brocker->getInterfaceImpl(type);
    }
    else {
        QMessageBox::warning(this, tr("error"), tr("can't cast brocker to type"));
        ILog::report("can't cast brocker to type");
        return false;
    }

    void* ptr;
    switch(type)
    {
    case IBrocker::PROBLEM:
        m_problem = reinterpret_cast<IProblem *>(interfaceImpl);
        ptr = (void*)m_problem;
        break;

    case IBrocker::SOLVER:
        m_solver = reinterpret_cast<ISolver *>(interfaceImpl);
        ptr = (void*)m_solver;
        break;

    default:
        QMessageBox::warning(this, tr("error"), tr("wrong type argument"));
        ILog::report("wrong type argument");
        return false;
    }

    if (!ptr) {
        QMessageBox::warning(this, tr("error"), tr("can't get implementation"));
        ILog::report("can't get implementation");
        return false;
    }

    //QMessageBox::warning(this, tr("info"), tr("load successfully"));
    return true;
}

void Controller::on_DrawButton_clicked()
{

}

void Controller::on_action_Insert_Row_triggered()
{
    insertRow();
}

void Controller::on_action_Delete_Row_triggered()
{
    deleteRow();
}

void Controller::on_pushButton_clicked()
{
    insertRow();
}

void Controller::on_pushButton_2_clicked()
{
    deleteRow();
}

void Controller::on_browseButton_clicked() //problem path
{
    QString fileName = QFileDialog::getOpenFileName(this, "Select 'problem' file", "", "DLL (*.dll)");
    ui->problemEdit->setText(fileName);
}

void Controller::closeEvent(QCloseEvent *event)
{
    int count = windows->count();
    for (int i = 0; i < count; i++){
        windows->at(i)->close();
    }
    this->close();
    event->accept();
}

void Controller::on_Controller_destroyed()
{

}

void Controller::on_SolveButton_released()
{
    QString problem_path = ui->problemEdit->text();
    if (problem_path.isEmpty()) {
        QMessageBox::warning(this, tr("error"),tr("problem path is empty"));
        ILog::report("problem path is empty");
        return;
    }

    if (!loadBrocker(IBrocker::PROBLEM, problem_path)) {
        ILog::report("Error loading problem");
        return;
    }

    if (!m_problem) {
        QMessageBox::warning(this, tr("error"),tr("No problem selected!"));
        ILog::report("No problem selected!");
        return;
    }
    QItemSelectionModel *selectionModel = ui->table->selectionModel();
    if (!selectionModel->hasSelection()) {
        QMessageBox::warning(this, tr("error"),tr("No solver selected!"));
        ILog::report("No solver selected!");
        return;
    }
    const int pathColumn = 2;
    QModelIndexList indexList = selectionModel->selectedRows(pathColumn);
    QVariant data = indexList.front().data();
    QString solver_path = data.toString();

    if (!loadBrocker(IBrocker::SOLVER, solver_path)) {
        ILog::report("Error loading solver");
        return;
    }
    ParamsDialog* dialog = new ParamsDialog(this, m_problem, m_solver);
    dialog->show();
    /*if (dialog.exec() != QDialog::Accepted) {
        QMessageBox::warning(this, tr("error"),tr("error opening dialog"));
        ILog::report("error opening dialog");
        return;
    }*/
}
