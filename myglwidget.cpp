#include <QtWidgets>
#include <QtOpenGL>
#include <GL/glu.h>
#include "myglwidget.h"
#include "help.hpp"
#include "chebyshev.hpp"
#include "auto_parallelization.hpp"
#define ARGC 6 //старое 5
MyGLWidget::MyGLWidget(QWidget *parent)
    : QGLWidget(QGLFormat(QGL::SampleBuffers), parent){
    xRot = 0;
    yRot = 0;
    zRot = 0;
}

MyGLWidget::~MyGLWidget(){
    free(F);
    free(cx);
    free(cy);
    free(Fx);
    free(Fy);
    free(TT);
    free(T);
}

int MyGLWidget::parse_command_line(){
    FILE *ff;
    QStringList args=QApplication::arguments();
    if(args.size()!= ARGC)
        return -3;
    ff=fopen(args.at(1).toLocal8Bit().constData(), "r");
    if(!ff)
        return -1;
    fscanf(ff, "%lf", &a);
    fscanf(ff, "%lf", &c);
    fscanf(ff, "%lf", &b);
    fscanf(ff, "%lf", &d);
    if(fabs(a-b)<1e-6 || fabs(c-d)<1e-6)
        return -2;
    fclose(ff);
    bool ok=true;
    nx=args.at(2).toInt(&ok);
    if(!ok)
        return -3;
    if(nx<3)
        return -2;
    ny=args.at(3).toInt(&ok);
    if(!ok)
        return -3;
    if(ny<3)
        return -2;
    k=args.at(4).toInt(&ok);
    if(!ok)
        return -3;
    if(k<0 || k>7)
        return -2;
    nthreads=args.at(5).toInt(&ok);
    if(!ok)
        return -3;
    if(nthreads < 1)
        return -2;
    allocate();
    chebyshevpoints(cx, cy, nx, ny, a, b, c, d);
    if (nthreads > 1){
        fill_Fx_Fy(Fx, cx, nx, a, b, Fy, cy, ny, c, d);
    }
    else {
        fill_Fx(Fx, cx, nx, a, b);
        fill_Fy(Fy, cy, ny, c, d);
    }
    change_func();
    interpolation_tensor(T, Fx, F, Fy, nx, ny, nthreads);
    Fill_TT(Fx, nx, T,Fy,ny, TT, nthreads);
    extrema_hunt();
    print_console();
    return 0;
}

void MyGLWidget::allocate(){
    if(F)
        free(F);
    if(cx)
        free(cx);
    if(cy)
        free(cy);
    if(Fx)
        free(Fx);
    if(Fy)
        free(Fy);
    if(T)
        free(T);
    if(TT)
        free(TT);
    F=(double*)malloc((nx+2)*(ny+2)*sizeof(double));
    cx=(double*)malloc((nx+2)*sizeof(double));
    cy=(double*)malloc((ny+2)*sizeof(double));
    Fx=(double*)malloc(nx*nx*sizeof(double));
    Fy=(double*)malloc(ny*ny*sizeof(double));
    T=(double*)malloc(nx*ny*sizeof(double));
    TT=(double*)malloc((nx+2)*(ny+2)*sizeof(double));
}

void MyGLWidget::print_console(){
    printf("format: %d\n", view_id);
    printf("area: [%lf;%lf]x[%lf;%lf]\n", a,b,c,d);
    printf("squeeze-stretch: %d\n", s);
    printf("points: %d %d\n", nx,ny);
    printf("disturbance: %d\n", p);
    printf("function absmax: %10.3e\n", absmax);
    printf("factic absmax: %10.3e\n", max(fabs(extr[0]),fabs(extr[1])));
    printf("extrema: %10.3e;%10.3e\n", extr[0], extr[1]);
    printf("angle around Oz: %d\n\n", zRot);
}

void MyGLWidget::change_func(){
    switch(k){
        case 0:
            f_name="k=0 f(x,y)=1";
            f=f0;
            break;
        case 1:
            f_name="k=1 f(x,y)=x";
            f=f1;
            break;
        case 2:
            f_name="k=2 f(x,y)=y";
            f=f2;
            break;
        case 3:
            f_name="k=3 f(x,y)=x+y";
            f=f3;
            break;
        case 4:
            f_name="k=4 f(x,y)=sqrt(x^2+y^2)";
            f=f4;
            break;
        case 5:
            f_name="k=5 f(x,y)=x^2+y^2";
            f=f5;
            break;
        case 6:
            f_name="k=6 f(x,y)=exp(x^2-y^2)";
            f=f6;
            break;
        case 7:
            f_name="k=7 f(x,y)=1/(25*(x^2+y^2)+1)";
            f=f7;
            break;
    }
    Fill_F(F, nx, ny, cx, cy, f, nthreads);
    max_z=max_matr(F,nx,ny);
    min_z=min_matr(F,nx,ny);
    absmax=max(fabs(max_z),fabs(min_z));
    if(p){
        F[(ny+2)*(nx/2+1)+ny/2+1]+=(0.1*absmax*p);
        if(F[(ny+2)*(nx/2+1)+ny/2+1]>max_z)
            max_z=F[(ny+2)*(nx/2+1)+ny/2+1];
        if(F[(ny+2)*(nx/2+1)+ny/2+1]<min_z)
            min_z=F[(ny+2)*(nx/2+1)+ny/2+1];
    }
}

void MyGLWidget::extrema_hunt(){
    if(view_id==0){
        extr[1]=max_z;
        extr[0]=min_z;
    }
    if(view_id==1){
        extr[1]=max_matr(TT,nx,ny);
        extr[0]=min_matr(TT,nx,ny);
    }
    if(view_id==2){
        extr[0]=F[0]-TT[0];
        extr[1]=F[0]-TT[0];
        for(int i=0; i<=nx+1; ++i){
            for(int j=0; j<=ny+1; ++j){
                if(F[i*(ny+2)+j]-TT[i*(ny+2)+j]>extr[1])
                    extr[1]=F[i*(ny+2)+j]-TT[i*(ny+2)+j];
                if(F[i*(ny+2)+j]-TT[i*(ny+2)+j]<extr[0])
                    extr[0]=F[i*(ny+2)+j]-TT[i*(ny+2)+j];
            }
        }
    }
}

void MyGLWidget::draw_area(){
    glColor3f(0.0,0.0,0.0);
    for(int i=0; i<=nx+1; ++i){
        glBegin(GL_LINE_STRIP);
        for(int j=0; j<=ny+1; ++j)
            glVertex3f(cx[i],cy[j],0);
        glEnd();
    }
    for(int j=0; j<=ny+1; ++j){
        glBegin(GL_LINE_STRIP);
        for(int i=0; i<=nx+1; ++i)
            glVertex3f(cx[i],cy[j],0);
        glEnd();
    }
}

void MyGLWidget::func_graph(){
    glColor3f(1.0,1.0,0.0);
    for(int i=0; i<=nx+1; ++i){
        glBegin(GL_LINE_STRIP);
        for(int j=0; j<=ny+1; ++j)
            glVertex3f(cx[i],cy[j],F[i*(ny+2)+j]);
        glEnd();
    }
    for(int j=0; j<=ny+1; ++j){
        glBegin(GL_LINE_STRIP);
        for(int i=0; i<=nx+1; ++i)
            glVertex3f(cx[i],cy[j],F[i*(ny+2)+j]);
        glEnd();
    }
}

void MyGLWidget::appr_graph(){
    glColor3f(1.0,0.0,1.0);
    for(int i=0; i<=nx+1; ++i){
        glBegin(GL_LINE_STRIP);
        for(int j=0; j<=ny+1; ++j)
            glVertex3f(cx[i],cy[j],TT[i*(ny+2)+j]);
        glEnd();
    }
    for(int j=0; j<=ny+1; ++j){
        glBegin(GL_LINE_STRIP);
        for(int i=0; i<=nx+1; ++i)
            glVertex3f(cx[i],cy[j],TT[i*(ny+2)+j]);
        glEnd();
    }
}

void MyGLWidget::err_graph(){
    glColor3f(0.0,1.0,1.0);
    for(int i=0; i<=nx+1; ++i){
        glBegin(GL_LINE_STRIP);
        for(int j=0; j<=ny+1; ++j)
            glVertex3f(cx[i],cy[j],F[i*(ny+2)+j]-TT[i*(ny+2)+j]);
        glEnd();
    }
    for(int j=0; j<=ny+1; ++j){
        glBegin(GL_LINE_STRIP);
        for(int i=0; i<=nx+1; ++i)
            glVertex3f(cx[i],cy[j],F[i*(ny+2)+j]-TT[i*(ny+2)+j]);
        glEnd();
    }
}

void MyGLWidget::press0(){
    change_func();
    interpolation_tensor(T, Fx, F, Fy, nx, ny, nthreads);
    Fill_TT(Fx, nx, T,Fy,ny, TT, nthreads);
}

void MyGLWidget::press23(){
    chebyshevpoints(cx, cy, nx, ny, a, b, c, d);
    if (nthreads > 1){
        fill_Fx_Fy(Fx, cx, nx, a, b, Fy, cy, ny, c, d);
    }
    else {
        fill_Fx(Fx, cx, nx, a, b);
        fill_Fy(Fy, cy, ny, c, d);
    }
    Fill_F(F, nx, ny, cx, cy, f, nthreads);
    max_z=max_matr(F,nx,ny);
    min_z=min_matr(F,nx,ny);
    absmax=max(fabs(max_z),fabs(min_z));
    if(p){
        F[(ny+2)*(nx/2+1)+ny/2+1]+=(0.1*absmax*p);
        if(F[(ny+2)*(nx/2+1)+ny/2+1]>max_z)
            max_z=F[(ny+2)*(nx/2+1)+ny/2+1];
        if(F[(ny+2)*(nx/2+1)+ny/2+1]<min_z)
            min_z=F[(ny+2)*(nx/2+1)+ny/2+1];
    }
    interpolation_tensor(T, Fx, F, Fy, nx, ny, nthreads);
    Fill_TT(Fx, nx, T,Fy,ny, TT, nthreads);
}

void MyGLWidget::press45(){
    allocate();
    chebyshevpoints(cx, cy, nx, ny, a, b, c, d);
    if (nthreads > 1){
        fill_Fx_Fy(Fx, cx, nx, a, b, Fy, cy, ny, c, d);
    }
    else {
        fill_Fx(Fx, cx, nx, a, b);
        fill_Fy(Fy, cy, ny, c, d);
    }
    Fill_F(F, nx, ny, cx, cy, f, nthreads);
    if(p){
        F[(ny+2)*(nx/2+1)+ny/2+1]+=(0.1*absmax*p);
        if(F[(ny+2)*(nx/2+1)+ny/2+1]>max_z)
            max_z=F[(ny+2)*(nx/2+1)+ny/2+1];
        if(F[(ny+2)*(nx/2+1)+ny/2+1]<min_z)
            min_z=F[(ny+2)*(nx/2+1)+ny/2+1];
    }
    interpolation_tensor(T, Fx, F, Fy, nx, ny, nthreads);
    Fill_TT(Fx, nx, T,Fy,ny, TT, nthreads);
}

void MyGLWidget::press67(){
    if(F[(ny+2)*(nx/2+1)+ny/2+1]>max_z)
        max_z=F[(ny+2)*(nx/2+1)+ny/2+1];
    if(F[(ny+2)*(nx/2+1)+ny/2+1]<min_z)
        min_z=F[(ny+2)*(nx/2+1)+ny/2+1];
    interpolation_tensor(T, Fx, F, Fy, nx, ny, nthreads);
    Fill_TT(Fx, nx, T,Fy,ny, TT, nthreads);
}

void MyGLWidget::printwindow(){
    //TODO:
    // switch to:   renderText(x,y,"....")
    qglColor(Qt::black);
    renderText(0, 15, f_name);
    renderText(10, 30, QString("format: %1").arg(view_id));
    renderText(10, 45, QString("scale: %1 [%2;%3]x[%4;%5]").arg(s).arg(a).arg(b).arg(c).arg(d));
    renderText(10, 60, QString("points: %1,%2").arg(nx).arg(ny));
    renderText(10, 75, QString("p: %1").arg(p));
    //if(k==0 && p==0 && view_id==0)
    //    renderText(10, 90, QString("absmax(fact): %1( %2)").arg(absmax).arg(absmax));
    //else
    renderText(10, 90, QString("absmax(fact): %1( %2)").arg(absmax).arg(max(fabs(extr[0]), fabs(extr[1]))));
    renderText(10, 105, QString("zRot: %1").arg(zRot));
}

QSize MyGLWidget::minimumSizeHint() const{
    return QSize(50, 50);
}

QSize MyGLWidget::sizeHint() const{
    return QSize(400, 400);
}

static void qNormalizeAngle(int &angle){
    while (angle < 0)
        angle += 360 * 16;
    while (angle > 360)
        angle -= 360 * 16;
}

void MyGLWidget::setXRotation(int angle){
    qNormalizeAngle(angle);
    if (angle != xRot) {
        xRot = angle;
        angle = angle>180 ? angle-360 : angle;
        emit xRotationChanged(angle);
        updateGL();
    }
}

void MyGLWidget::setYRotation(int angle){
    qNormalizeAngle(angle);
    if (angle != yRot) {
        yRot = angle;
        angle = angle>180 ? angle-360 : angle;
        emit yRotationChanged(angle);
        updateGL();
    }
}

void MyGLWidget::setZRotation(int angle){
    qNormalizeAngle(angle);
    if (angle != zRot) {
        zRot = angle;
        angle = angle>180 ? angle-360 : angle;
        emit zRotationChanged(angle);
        updateGL();
    }
}

void MyGLWidget::initializeGL(){
    qglClearColor(Qt::white);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);
    static GLfloat lightPosition[4] = { 0, 0, 10, 1.0 };
    glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
    glEnable(GL_LIGHT0);
}

void MyGLWidget::setProjection() {
    // вычислим сдвиг и масштаб так, чтобы куб, содержащий график
    // вписывался в окно [-1, 1]x[-1, 1]
    double modelM[16],projM[16];
    int view[4]={-1,-1,2,2};
    double verts[8][3];
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glGetDoublev(GL_MODELVIEW_MATRIX,modelM);
    glGetDoublev(GL_PROJECTION_MATRIX,projM);
    //glGetIntegerv(GL_VIEWPORT,view);
    // ищем проекции вершин куба [a,b]x[c,d]x[minz,maxz] на плоскость экрана
    /*gluProject(a,c,extr[0],modelM,projM,view,
            &verts[0][0],&verts[0][1],&verts[0][2]);
    gluProject(b,c,extr[0],modelM,projM,view,
            &verts[1][0],&verts[1][1],&verts[1][2]);
    gluProject(a,d,extr[0],modelM,projM,view,
            &verts[2][0],&verts[2][1],&verts[2][2]);
    gluProject(b,d,extr[0],modelM,projM,view,
            &verts[3][0],&verts[3][1],&verts[3][2]);
    gluProject(a,c,extr[1],modelM,projM,view,
            &verts[4][0],&verts[4][1],&verts[4][2]);
    gluProject(b,c,extr[1],modelM,projM,view,
            &verts[5][0],&verts[5][1],&verts[5][2]);
    gluProject(a,d,extr[1],modelM,projM,view,
            &verts[6][0],&verts[6][1],&verts[6][2]);
    gluProject(b,d,extr[1],modelM,projM,view,
            &verts[7][0],&verts[7][1],&verts[7][2]);*/
    if(view_id==2){
        gluProject(a,c,extr[0],modelM,projM,view,
            &verts[0][0],&verts[0][1],&verts[0][2]);
        gluProject(b,c,extr[0],modelM,projM,view,
            &verts[1][0],&verts[1][1],&verts[1][2]);
        gluProject(a,d,extr[0],modelM,projM,view,
            &verts[2][0],&verts[2][1],&verts[2][2]);
        gluProject(b,d,extr[0],modelM,projM,view,
            &verts[3][0],&verts[3][1],&verts[3][2]);
        gluProject(a,c,extr[1],modelM,projM,view,
            &verts[4][0],&verts[4][1],&verts[4][2]);
        gluProject(b,c,extr[1],modelM,projM,view,
            &verts[5][0],&verts[5][1],&verts[5][2]);
        gluProject(a,d,extr[1],modelM,projM,view,
            &verts[6][0],&verts[6][1],&verts[6][2]);
        gluProject(b,d,extr[1],modelM,projM,view,
            &verts[7][0],&verts[7][1],&verts[7][2]);
    }else{
        gluProject(a,c,min(extr[0],0),modelM,projM,view,
            &verts[0][0],&verts[0][1],&verts[0][2]);
        gluProject(b,c,min(extr[0],0),modelM,projM,view,
            &verts[1][0],&verts[1][1],&verts[1][2]);
        gluProject(a,d,min(extr[0],0),modelM,projM,view,
            &verts[2][0],&verts[2][1],&verts[2][2]);
        gluProject(b,d,min(extr[0],0),modelM,projM,view,
            &verts[3][0],&verts[3][1],&verts[3][2]);
        gluProject(a,c,max(extr[1],0),modelM,projM,view,
            &verts[4][0],&verts[4][1],&verts[4][2]);
        gluProject(b,c,max(extr[1],0),modelM,projM,view,
            &verts[5][0],&verts[5][1],&verts[5][2]);
        gluProject(a,d,max(extr[1],0),modelM,projM,view,
            &verts[6][0],&verts[6][1],&verts[6][2]);
        gluProject(b,d,max(extr[1],0),modelM,projM,view,
            &verts[7][0],&verts[7][1],&verts[7][2]);
    }
    // ищем квадрат в плоскости экрана куда они вписываются
    double minX=verts[0][0], maxX=verts[0][0],
            minY=verts[0][1], maxY=verts[0][1],
            minZ=verts[0][2], maxZ=verts[0][2];
    for (int i=1; i<8; i++) {
        if (verts[i][0]<minX) minX=verts[i][0];
        if (verts[i][1]<minY) minY=verts[i][1];
        if (verts[i][2]<minZ) minZ=verts[i][2];
        if (verts[i][0]>maxX) maxX=verts[i][0];
        if (verts[i][1]>maxY) maxY=verts[i][1];
        if (verts[i][2]>maxZ) maxZ=verts[i][2];
    }
    double sz=maxX-minX;
    if (maxY-minY>sz) sz=maxY-minY;
    sz*=1.1;
    // "центр масс"
    double sx=(maxX+minX)/2,sy=(maxY+minY)/2;
    glOrtho(sx-sz/2, sx+sz/2, sy-sz/2, sy+sz/2, -20, +20);
    glMatrixMode(GL_MODELVIEW);
}

void MyGLWidget::paintGL(){
    printwindow();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    // TODO:   [a,b]x[c,d]   Scale + Translate --> [-1,1]x[-1,1]
    glTranslatef(0.0, 0.0, -10.0);
    glRotatef(zRot , 0.0, 0.0, 1.0);
    glRotatef(xRot , 1.0, 0.0, 0.0);
    glRotatef(yRot , 0.0, 1.0, 0.0);
    setProjection();
    glLineWidth(2.0);
    glBegin(GL_LINES);
    glColor3d(1.0,0.0,0);
    glVertex3d(a,0,0);
    glVertex3d(b,0,0);
    glColor3d(0.0,1.0,0);
    glVertex3d(0,c,0);
    glVertex3d(0,d,0);
    glColor3d(0,0.0,1.0);
    glVertex3d(0,0,min(0,extr[0]));
    glVertex3d(0,0,max(0,extr[1]));
    glEnd();
    if(view_id!=2)
        draw_area();
    if(view_id==0)
        func_graph();
    if(view_id==1)
        appr_graph();
    if(view_id==2)
        err_graph();
    printwindow();
}

void MyGLWidget::resizeGL(int width, int height){
    glViewport(0, 0, width, height);
    /*double aspect=height*1.0/width;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    // TODO: if aspect>1 ....
    glOrtho(-1, +1, -1*aspect, +1*aspect, 1.0, 15.0);
    glMatrixMode(GL_MODELVIEW);*/
}

void MyGLWidget::mousePressEvent(QMouseEvent *event){
    lastPos = event->pos();
}

void MyGLWidget::mouseMoveEvent(QMouseEvent *event){
    int dx = event->x() - lastPos.x();
    int dy = event->y() - lastPos.y();
    if (event->buttons() & Qt::LeftButton) {
        setXRotation(xRot + dy);
        setYRotation(yRot + dx);
    } else if (event->buttons() & Qt::RightButton) {
        setXRotation(xRot + dy);
        setZRotation(zRot + dx);
    }
    lastPos = event->pos();
}

void MyGLWidget::keyPressEvent(QKeyEvent* e){
    switch (e->key()){
        case Qt::Key_0:
            k=(k+1)%8;
            press0();
            break;
        case Qt::Key_1:
            view_id=(view_id+1)%3;
            break;
        case Qt::Key_2:
            { double w=(b-a)/2,h=(d-c)/2;
            a=a-w;
            b=b+w;
            c=c-h;
            d=d+h;}
            ++s;
            press23();
            break;
        case Qt::Key_3:
            { double w=(b-a)/4,h=(d-c)/4;
            a=a+w;
            b=b-w;
            c=c+h;
            d=d-h;}
            --s;
            press23();
            break;
        case Qt::Key_4:
            nx*=2;
            ny*=2;
            press45();
            break;
        case Qt::Key_5:
            nx/=2;
            ny/=2;
            press45();
            break;
        case Qt::Key_6:
            ++p;
            F[(ny+2)*(nx/2+1)+ny/2+1]+=(0.1*absmax);
            press67();
            break;
        case Qt::Key_7:
            --p;
            F[(ny+2)*(nx/2+1)+ny/2+1]-=(0.1*absmax);
            press67();
            break;
        case Qt::Key_8:
            setZRotation(zRot+15);
            break;
        case Qt::Key_9:
            setZRotation(zRot-15);
            break;
    }
    if (e->key()==Qt::Key_X &&
            (e->modifiers() & Qt::ControlModifier)){
        close();
    }
    e->accept();
    extrema_hunt();
    print_console();
    updateGL();
}
