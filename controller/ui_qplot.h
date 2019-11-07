/********************************************************************************
** Form generated from reading UI file 'qplot.ui'
**
** Created: Wed 4. Sep 23:36:34 2019
**      by: Qt User Interface Compiler version 4.7.4
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_QPLOT_H
#define UI_QPLOT_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QDoubleSpinBox>
#include <QtGui/QGridLayout>
#include <QtGui/QGroupBox>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QPushButton>
#include <QtGui/QSpinBox>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>
#include "qcustomplot.h"

QT_BEGIN_NAMESPACE

class Ui_QPlot
{
public:
    QGroupBox *qraphic;
    QVBoxLayout *verticalLayout_4;
    QGridLayout *gridLayout;
    QLabel *labelAxis;
    QLabel *labelLB;
    QLabel *labelRB;
    QSpinBox *spinAxis;
    QDoubleSpinBox *spinLeftBorder;
    QDoubleSpinBox *spinRightBorder;
    QPushButton *drawButton;
    QCustomPlot *plot;

    void setupUi(QWidget *QPlot)
    {
        if (QPlot->objectName().isEmpty())
            QPlot->setObjectName(QString::fromUtf8("QPlot"));
        QPlot->resize(714, 467);
        qraphic = new QGroupBox(QPlot);
        qraphic->setObjectName(QString::fromUtf8("qraphic"));
        qraphic->setGeometry(QRect(10, 10, 701, 441));
        QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(qraphic->sizePolicy().hasHeightForWidth());
        qraphic->setSizePolicy(sizePolicy);
        verticalLayout_4 = new QVBoxLayout(qraphic);
        verticalLayout_4->setObjectName(QString::fromUtf8("verticalLayout_4"));
        gridLayout = new QGridLayout();
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        labelAxis = new QLabel(qraphic);
        labelAxis->setObjectName(QString::fromUtf8("labelAxis"));

        gridLayout->addWidget(labelAxis, 0, 0, 1, 1);

        labelLB = new QLabel(qraphic);
        labelLB->setObjectName(QString::fromUtf8("labelLB"));

        gridLayout->addWidget(labelLB, 1, 0, 1, 1);

        labelRB = new QLabel(qraphic);
        labelRB->setObjectName(QString::fromUtf8("labelRB"));

        gridLayout->addWidget(labelRB, 2, 0, 1, 1);

        spinAxis = new QSpinBox(qraphic);
        spinAxis->setObjectName(QString::fromUtf8("spinAxis"));
        spinAxis->setMaximum(0);

        gridLayout->addWidget(spinAxis, 0, 1, 1, 1);

        spinLeftBorder = new QDoubleSpinBox(qraphic);
        spinLeftBorder->setObjectName(QString::fromUtf8("spinLeftBorder"));
        spinLeftBorder->setMinimum(-99.99);

        gridLayout->addWidget(spinLeftBorder, 1, 1, 1, 1);

        spinRightBorder = new QDoubleSpinBox(qraphic);
        spinRightBorder->setObjectName(QString::fromUtf8("spinRightBorder"));
        spinRightBorder->setMinimum(-99.99);
        spinRightBorder->setValue(5);

        gridLayout->addWidget(spinRightBorder, 2, 1, 1, 1);


        verticalLayout_4->addLayout(gridLayout);

        drawButton = new QPushButton(qraphic);
        drawButton->setObjectName(QString::fromUtf8("drawButton"));

        verticalLayout_4->addWidget(drawButton);

        plot = new QCustomPlot(qraphic);
        plot->setObjectName(QString::fromUtf8("plot"));

        verticalLayout_4->addWidget(plot);


        retranslateUi(QPlot);

        QMetaObject::connectSlotsByName(QPlot);
    } // setupUi

    void retranslateUi(QWidget *QPlot)
    {
        QPlot->setWindowTitle(QApplication::translate("QPlot", "Form", 0, QApplication::UnicodeUTF8));
        qraphic->setTitle(QApplication::translate("QPlot", "Graphic", 0, QApplication::UnicodeUTF8));
        labelAxis->setText(QApplication::translate("QPlot", "Axis", 0, QApplication::UnicodeUTF8));
        labelLB->setText(QApplication::translate("QPlot", "Left border", 0, QApplication::UnicodeUTF8));
        labelRB->setText(QApplication::translate("QPlot", "Right border", 0, QApplication::UnicodeUTF8));
        drawButton->setText(QApplication::translate("QPlot", "Draw", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class QPlot: public Ui_QPlot {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_QPLOT_H
