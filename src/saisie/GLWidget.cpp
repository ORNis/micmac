#include "GLWidget.h"

//Min and max zoom ratio (relative)
const float GL_MAX_ZOOM_RATIO = 1.0e6f;
const float GL_MIN_ZOOM_RATIO = 1.0e-6f;

//invalid GL list index
const GLuint GL_INVALID_LIST_ID = (~0);

using namespace Cloud_;
using namespace std;

GLWidget::GLWidget(QWidget *parent, cData *data) : QGLWidget(parent)
  , m_rw(1.f)
  , m_rh(1.f)
  , m_font(font())
  , m_bDrawAxis(false)
  , m_bDrawBall(true)
  , m_bDrawCams(true)
  , m_bDrawMessages(true)
  , m_bDrawBbox(false)
  , m_bObjectCenteredView(true)
  , m_bPolyIsClosed(false)
  , m_interactionMode(TRANSFORM_CAMERA)
  , m_bFirstAction(true)
  , m_trihedronGLList(GL_INVALID_LIST_ID)
  , m_ballGLList(GL_INVALID_LIST_ID)
  , m_textureImage(GL_INVALID_LIST_ID)
  , m_textureMask(GL_INVALID_LIST_ID)
  , m_nbGLLists(0)
  , m_params(ViewportParameters())
  , m_Data(data)
  , m_speed(2.5f)
  , m_bDisplayMode2D(false)
  , m_Click(0)
  , m_sqr_radius(2500)
  , m_vertexbuffer(QGLBuffer::VertexBuffer)
  , _frameCount(0)
  , _previousTime(0)
  , _currentTime(0)
  , _idx(-1)
  , _fps(0.0f)
  , _m_g_mouseLeftDown(false)
  , _m_g_mouseMiddleDown(false)
  , _m_g_mouseRightDown(false)
  , _mask(NULL)
{
    _m_g_rotationMatrix[0] = _m_g_rotationMatrix[4] = _m_g_rotationMatrix[8] = 1;
    _m_g_rotationMatrix[1] = _m_g_rotationMatrix[2] = _m_g_rotationMatrix[3] = 0;
    _m_g_rotationMatrix[5] = _m_g_rotationMatrix[6] = _m_g_rotationMatrix[7] = 0;  

    _time.start();

    setFocusPolicy(Qt::StrongFocus);

    //drag & drop handling
    setAcceptDrops(true);

    m_glPosition[0] = m_glPosition[1] = 0.f;

    _mvmatrix   = new GLdouble[16];
    _projmatrix = new GLdouble[16];
    _glViewport = new GLint[4];

    m_font.setPointSize(10);
}

GLWidget::~GLWidget()
{
    if (m_trihedronGLList != GL_INVALID_LIST_ID)
    {
        glDeleteLists(m_trihedronGLList,1);
        m_trihedronGLList = GL_INVALID_LIST_ID;
    }
    if (m_ballGLList != GL_INVALID_LIST_ID)
    {
        glDeleteLists(m_ballGLList,1);
        m_ballGLList = GL_INVALID_LIST_ID;
    }
    if (m_textureImage != GL_INVALID_LIST_ID)
    {
        glDeleteLists(m_textureImage,1);
        m_textureImage = GL_INVALID_LIST_ID;
    }
}

bool GLWidget::eventFilter(QObject* object,QEvent* event)
{
    if(event->type() == QEvent::MouseMove)
    {
        QMouseEvent *mouseEvent = static_cast<QMouseEvent*>(event);
        m_lastPos = mouseEvent->pos();

        update();

        if (m_interactionMode == SELECTION)
        {
            if(!m_bPolyIsClosed)
            {
                int sz = m_polygon.size();

                if (m_Data->NbImages())
                {
                    if (sz == 0)
                        return false;
                    else if (sz == 1)
                        m_polygon.push_back(WindowToImage(m_lastPos));
                    else
                        //replace last point by the current one
                        m_polygon[sz-1] = WindowToImage(m_lastPos);

                }
                else
                {
                    if (sz == 0)
                        return false;
                    else if (sz == 1)
                        m_polygon.push_back(m_lastPos);
                    else
                        //replace last point by the current one
                        m_polygon[sz-1] = m_lastPos;
                }
            }
            else
            {
                if(m_polygon.size())
                {
                    QKeyEvent *keyEvent = static_cast<QKeyEvent*>(event);
                    if (keyEvent->modifiers().testFlag(Qt::ShiftModifier))
                    {
                        fillPolygon2();
                    }
                    else
                    {
                        if (m_Click == 1)
                        {
                            m_polygon2.clear();
                            if ((_idx >0 ) && (_idx < m_polygon.size()-1))
                            {
                                m_polygon2.push_back(m_polygon[_idx-1]);
                                m_polygon2.push_back(m_lastPos);
                                m_polygon2.push_back(m_polygon[_idx+1]);
                            }
                            else
                            {
                                if (_idx == 0)
                                {
                                    m_polygon2.push_back(m_polygon[m_polygon.size()-1]);
                                    m_polygon2.push_back(m_lastPos);
                                    m_polygon2.push_back(m_polygon[1]);
                                }
                                if (_idx == m_polygon.size()-1)
                                {
                                    m_polygon2.push_back(m_polygon[m_polygon.size()-2]);
                                    m_polygon2.push_back(m_lastPos);
                                    m_polygon2.push_back(m_polygon[0]);
                                }
                            }
                        }
                        else
                        {
                            findClosestPoint();
                        }
                    }
                }
            }

            update();

        }

        return true;
    }  
    else
    {
          return QObject::eventFilter(object,event);
    }
}

void GLWidget::initializeGL()
{
    if (m_bGLInitialized)
        return;

    glEnable( GL_DEPTH_TEST );

    m_bGLInitialized = true;
}

void GLWidget::resizeGL(int width, int height)
{
    if (width==0 || height==0) return;

    m_glRatio  = (float) width/height;

    glViewport( 0, 0, width, height );
    glGetIntegerv (GL_VIEWPORT, _glViewport);

    if (m_Data->NbImages())
        zoomFit();
}

//-------------------------------------------------------------------------
// Computes the frames rate
//-------------------------------------------------------------------------
void GLWidget::computeFPS()
{
    //  Increase frame count
    _frameCount++;

    _currentTime = _time.elapsed();

    //  Compute elapsed time
    int deltaTime = _currentTime - _previousTime;

    if(deltaTime > 1000)
    {
        //  compute the number of frames per second
        _fps = _frameCount / (deltaTime / 1000.0f);

        //  Set time
        _previousTime = _currentTime;

        //  Reset frame count
        _frameCount = 0;

        if (_fps > 1e-3)
        {
            m_messageFPS = "fps: " + QString::number(_fps);
        }
    }
}

//draw a unit circle in a given plane (0=YZ, 1 = XZ, 2=XY)
void glDrawUnitCircle(uchar dim, float cx, float cy, float r, int steps = 64)
{
    float theta = 2.f * PI / float(steps);
    float c = cosf(theta);//precalculate the sine and cosine
    float s = sinf(theta);
    float t;

    float x = r;//we start at angle = 0
    float y = 0;

    uchar dimX = (dim<2 ? dim+1 : 0);
    uchar dimY = (dimX<2 ? dimX+1 : 0);

    GLfloat P[3];

    for (int i=0;i<3;++i) P[i] = 0.0f;

    glBegin(GL_LINE_LOOP);
    for(int ii = 0; ii < steps; ii++)
    {
        P[dimX] = x + cx;
        P[dimY] = y + cy;
        glVertex3fv(P);

        //apply the rotation matrix
        t = x;
        x = c * x - s * y;
        y = s * t + c * y;
    }

    glEnd();
}

void GLWidget::drawQuad(GLfloat originX, GLfloat originY, GLfloat glh, GLfloat glw)
{
    glBegin(GL_QUADS);
    {
        glTexCoord2f(0.0f, 0.0f);
        glVertex2f(originX, originY);
        glTexCoord2f(1.0f, 0.0f);
        glVertex2f(originX+glw, originY);
        glTexCoord2f(1.0f, 1.0f);
        glVertex2f(originX+glw, originY+glh);
        glTexCoord2f(0.0f, 1.0f);
        glVertex2f(originX, originY+glh);
    }
    glEnd();
}

void GLWidget::drawQuad(GLfloat originX, GLfloat originY, GLfloat glh, GLfloat glw, QColor color)
{
    glColor4f(color.redF(),color.greenF(),color.blueF(),color.alphaF());
    drawQuad(originX,originY,glh,glw);
}

void GLWidget::drawQuad(GLfloat originX, GLfloat originY, GLfloat glh, GLfloat glw, GLuint idTexture)
{
    glEnable(GL_TEXTURE_2D);
    glBindTexture( GL_TEXTURE_2D, idTexture );
    drawQuad(originX,originY,glh,glw);
    glBindTexture( GL_TEXTURE_2D, 0);
    glDisable(GL_TEXTURE_2D);
}

void GLWidget::enableOptionLine()
{
    glEnable (GL_LINE_SMOOTH);
    glEnable (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glHint (GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
}

void GLWidget::disableOptionLine()
{
    glDisable(GL_BLEND);
    glDisable (GL_LINE_SMOOTH);
}

void GLWidget::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE,GL_ZERO);

    //gradient color background
    drawGradientBackground();
    //we clear background
    glClear(GL_DEPTH_BUFFER_BIT);

    glDisable(GL_BLEND);

    if (m_Data->NbImages())
    {

        glEnable(GL_BLEND);
        glBlendFunc(GL_ONE,GL_ZERO);

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();

        glDisable(GL_ALPHA_TEST);
        glDisable(GL_DEPTH_TEST);

        GLfloat glw = 2*m_rw;
        GLfloat glh = 2*m_rh;

        glPushMatrix();
        glMultMatrixd(_projmatrix);

        if(_projmatrix[0] != m_params.zoom)
        {
            GLint recal;
            GLdouble wx, wy, wz;

            recal = _glViewport[3] - (GLint) _m_lastPosZoom.y()- 1;

            gluUnProject ((GLdouble) _m_lastPosZoom.x(), (GLdouble) recal, 1.0,
                          _mvmatrix, _projmatrix, _glViewport, &wx, &wy, &wz);

            glTranslatef(wx,wy,0);
            glScalef(m_params.zoom/_projmatrix[0], m_params.zoom/_projmatrix[0], 1.0);
            glTranslatef(-wx,-wy,0);
        }

        glTranslatef(m_glPosition[0],m_glPosition[1],0);

        m_glPosition[0] = 0;
        m_glPosition[1] = 0;

        glGetDoublev (GL_PROJECTION_MATRIX, _projmatrix);
        drawQuad(0, 0, glh, glw,QColor(255,255,255));

        if(_mask != NULL && !_m_g_mouseMiddleDown)
        {

            drawQuad(0,0,glh,glw,m_textureMask );
            glBlendFunc(GL_ONE,GL_ONE);

            drawQuad(0, 0, glh, glw,QColor(128,128,128));
            glBlendFunc(GL_DST_COLOR,GL_SRC_COLOR);
        }

        drawQuad(0,0,glh,glw,m_textureImage );

        glPopMatrix();

        glDisable(GL_BLEND);
        glEnable(GL_DEPTH_TEST);
        glEnable(GL_ALPHA_TEST);
        glMatrixMode(GL_MODELVIEW);

        //Affichage du zoom et des coordonnées image
        if (m_bDrawMessages)
        {
            QPointF ptImg = WindowToImage(m_lastPos);

            glColor3f(1.f,1.f,1.f);

            renderText(10, _glViewport[3] - m_font.pointSize(), QString::number(m_params.zoom*100,'f',1) + "%", m_font);

            if  ((m_interactionMode == SELECTION)&&(ptImg.x()>=0)&&(ptImg.y()>=0)&&(ptImg.x()<_glImg.width())&&(ptImg.y()<_glImg.height()))
                renderText(_glViewport[2] - 120, _glViewport[3] - m_font.pointSize(), QString::number(ptImg.x(),'f',1) + ", " + QString::number(_glImg.height()-ptImg.y(),'f',1) + " px", m_font);

        }
    }
    else
    {


        zoom();

        static GLfloat trans44[16], rot44[16], tmp[16];
        m33_to_m44( _m_g_rotationMatrix, rot44 );
        setTranslate_m3(  m_params.m_translationMatrix, trans44 );

        //mult( trans44, rot44, tmp );
        mult( rot44, trans44, tmp );
        transpose( tmp, _m_g_glMatrix );
        glLoadMatrixf( _m_g_glMatrix );

        if (m_Data->NbClouds())
        {
            glEnable(GL_DEPTH_TEST);

            glEnableClientState(GL_VERTEX_ARRAY);
            glEnableClientState(GL_COLOR_ARRAY);

            m_vertexbuffer.bind();
            glVertexPointer(3, GL_FLOAT, 0, NULL);
            m_vertexbuffer.release();

            m_vertexColor.bind();
            glColorPointer(3, GL_FLOAT, 0, NULL);
            m_vertexColor.release();

            glDrawArrays( GL_POINTS, 0, m_Data->getCloud(0)->size()*3 );

            glDisableClientState(GL_VERTEX_ARRAY);
            glDisableClientState(GL_COLOR_ARRAY);

            glDisable(GL_DEPTH_TEST);
        }

        enableOptionLine();

        if (m_bDrawBall) drawBall();
        else if (m_bDrawAxis) drawAxis();

        if (m_bDrawCams) drawCams();

        if (m_bDrawBbox) drawBbox();

        disableOptionLine();

        if (m_Data->NbClouds()&& m_bDrawMessages)
        {
            computeFPS();

            glColor4f(0.8f,0.9f,1.0f,0.9f);

            renderText(10, _glViewport[3]- m_font.pointSize(), m_messageFPS,m_font);
        }
    }

    if (m_interactionMode == SELECTION)
    {
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(0,_glViewport[2],_glViewport[3],0,-1,1);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        glDisable(GL_DEPTH_TEST);

        enableOptionLine();

        drawPolygon();

        glLineWidth(m_params.LineWidth);
        disableOptionLine();

        drawPointAndSegments();

        glEnable(GL_DEPTH_TEST);

        glPopMatrix(); // restore modelview
    }

    //current messages (if valid)
    if (!m_messagesToDisplay.empty())
    {
        //Some versions of Qt seem to need glColorf instead of glColorub! (see https://bugreports.qt-project.org/browse/QTBUG-6217)
        glColor3f(1.f,1.f,1.f);

        int lc_currentHeight = _glViewport[3] - m_font.pointSize()*m_messagesToDisplay.size(); //lower center
        int uc_currentHeight = 10;            //upper center

        std::list<MessageToDisplay>::iterator it = m_messagesToDisplay.begin();
        while (it != m_messagesToDisplay.end())
        {
            switch(it->position)
            {
            case LOWER_LEFT_MESSAGE:
            {
                renderText(10, lc_currentHeight, it->message,m_font);
                int messageHeight = QFontMetrics(m_font).height();
                lc_currentHeight -= (messageHeight*5)/4; //add a 25% margin
            }
                break;
            case LOWER_CENTER_MESSAGE:
            {
                QRect rect = QFontMetrics(m_font).boundingRect(it->message);
                renderText((_glViewport[2]-rect.width())/2, lc_currentHeight, it->message,m_font);
                int messageHeight = QFontMetrics(m_font).height();
                lc_currentHeight += (messageHeight*5)/4; //add a 25% margin
            }
                break;
            case UPPER_CENTER_MESSAGE:
            {
                QRect rect = QFontMetrics(m_font).boundingRect(it->message);
                renderText((_glViewport[2]-rect.width())/2, uc_currentHeight+rect.height(), it->message,m_font);
                uc_currentHeight += (rect.height()*5)/4; //add a 25% margin
            }
                break;
            case SCREEN_CENTER_MESSAGE:
            {
                m_font.setPointSize(12);
                QRect rect = QFontMetrics(m_font).boundingRect(it->message);
                renderText((_glViewport[2]-rect.width())/2, (_glViewport[3]-rect.height())/2, it->message,m_font);
                m_font.setPointSize(10);
            }
            }

            ++it;
        }
    }
}

void GLWidget::mousePressEvent(QMouseEvent *event)
{
    m_lastPos = event->pos();

    if ( event->button() == Qt::LeftButton )
    {
        _m_g_mouseLeftDown = true;

        if (m_interactionMode == SELECTION)
        {
            if (m_Data->NbImages()||m_Data->NbClouds()||m_Data->NbCameras())
            {
                if (!m_bPolyIsClosed)
                {
                    if (m_Data->NbImages())
                    {
                          QPointF ptImg = WindowToImage(m_lastPos);

                          if (m_polygon.size() >= 1)
                               m_polygon[m_polygon.size()-1] = ptImg;

                          m_polygon.push_back(ptImg);
                    }
                    else
                    {
                        if (m_polygon.size() >= 1)
                            m_polygon[m_polygon.size()-1] = m_lastPos;

                        m_polygon.push_back(m_lastPos);
                    }
                }
                else
                {
                    if (event->modifiers().testFlag(Qt::ShiftModifier))
                    {
                        if ((m_polygon.size() >=2) && m_polygon2.size() && m_bPolyIsClosed)
                        {
                            // modify polygon...
                            int idx = -1;

                            for (int i=0;i<m_polygon.size();++i)
                            {
                                if (m_polygon[i] == m_polygon2[0]) idx = i;
                            }

                            if (idx >=0) m_polygon.insert(idx+1, m_polygon2[1]);
                        }

                        m_polygon2.clear();
                        update();
                    }
                    else if (_idx != -1)
                        m_Click++;
                }
            }
        }
    }
    else if (event->button() == Qt::RightButton)
    {
        if (m_interactionMode == TRANSFORM_CAMERA)
            _m_g_mouseRightDown = true;
        else if ((_idx >=0)&&(_idx<m_polygon.size()))
        {
            m_polygon.remove(_idx);

            findClosestPoint();

            if (m_polygon.size() < 2) m_bPolyIsClosed = false;

            update();
        }
        else
        {
            closePolyline();
            update();
        }
    }
    else if (event->button() == Qt::MiddleButton)
    {

        if (m_interactionMode == TRANSFORM_CAMERA)
            _m_g_mouseMiddleDown = true;

        _m_lastPosZoom = m_lastPos;
        update();
    }
}

void GLWidget::mouseReleaseEvent(QMouseEvent *event)
{
    if ( event->button() == Qt::LeftButton )
    {
        _m_g_mouseLeftDown = false;

        if ((m_Click >=1) &&(_idx>=0)&&m_polygon2.size())
        {
            m_polygon[_idx] = m_polygon2[1];

            m_polygon2.clear();
            m_Click = 0;

            update();
        }

        if ((m_Click >=1))
        {
            findClosestPoint();

            update();
        }

    }
    if ( event->button() == Qt::RightButton  )
    {
        _m_g_mouseRightDown = false;
    }
    if ( event->button() == Qt::MiddleButton  )
    {
        _m_g_mouseMiddleDown = false;
        update();
    }
}

void GLWidget::keyPressEvent(QKeyEvent* event)
{
    if(event->modifiers().testFlag(Qt::ControlModifier))
    {
        if(event->key() == Qt::Key_1)    zoomFactor(50);
        else if(event->key() == Qt::Key_2)    zoomFactor(25);
    }
    else
    {
        switch(event->key())
        {
        case Qt::Key_Escape:
            clearPolyline();
            break;
        case Qt::Key_1:
            zoomFactor(100);
            break;
        case Qt::Key_2:
            zoomFactor(200);
            break;
        case Qt::Key_4:
            zoomFactor(400);
            break;
        case Qt::Key_9:
            zoomFit();
            break;
        case Qt::Key_Plus:
            if (m_bDisplayMode2D)
                setZoom(m_params.zoom*1.5f);
            else
                ptSizeUp(true);
            break;
        case Qt::Key_Minus:
            if (m_bDisplayMode2D)
                setZoom(m_params.zoom/1.5f);
            else
                ptSizeUp(false);
            break;
        default:
            event->ignore();
            break;
        }
    }
    update();
}

void GLWidget::keyReleaseEvent(QKeyEvent* event)
{
    if  (event->key() == Qt::Key_Shift)
    {
        m_polygon2.clear();
    }

    if  (event->key() == Qt::Key_Control)
    {
        m_polygon2.clear();
        m_Click = 0;
    }

     findClosestPoint();
}

void GLWidget::setBufferGl(bool onlyColor)
{
    if(m_vertexbuffer.isCreated() && !onlyColor)
        m_vertexbuffer.destroy();
    if(m_vertexColor.isCreated())
        m_vertexColor.destroy();

    int sizeClouds = m_Data->getSizeClouds();

    if (sizeClouds == 0) return;

    GLfloat* vertices = NULL, *colors = NULL;

    if(!onlyColor)
        vertices   = new GLfloat[sizeClouds*3];

    colors     = new GLfloat[sizeClouds*3];
    int pitchV = 0;

    for (int aK=0; aK < m_Data->NbClouds();++aK){

        Cloud* pcloud = m_Data->getCloud(aK);
        uint sizeCloud = pcloud->size();

        for(int bK=0; bK< (int)sizeCloud; bK++)
        {
            Vertex vert = pcloud->getVertex(bK);
            QColor colo = vert.getColor();
            if(!onlyColor)
            {
                vertices[pitchV+bK*3 + 0 ] = vert.x();
                vertices[pitchV+bK*3 + 1 ] = vert.y();
                vertices[pitchV+bK*3 + 2 ] = vert.z();
            }
            if(vert.isVisible())
            {
                colors[pitchV+bK*3 + 0 ]   = colo.redF();
                colors[pitchV+bK*3 + 1 ]   = colo.greenF();
                colors[pitchV+bK*3 + 2 ]   = colo.blueF();
            }
            else
            {
                colors[pitchV+bK*3 + 0 ]   = colo.redF()   *0.7;
                colors[pitchV+bK*3 + 1 ]   = colo.greenF() *0.6;
                colors[pitchV+bK*3 + 2 ]   = colo.blueF()  *0.8;
            }
        }

        pitchV += sizeCloud;
    }

    if(!onlyColor)
    {
        m_vertexbuffer.create();
        m_vertexbuffer.setUsagePattern(QGLBuffer::StaticDraw);
        m_vertexbuffer.bind();
        m_vertexbuffer.allocate(vertices, sizeClouds* 3 * sizeof(GLfloat));
        m_vertexbuffer.release();
    }

    m_vertexColor.create();
    m_vertexColor.setUsagePattern(QGLBuffer::StaticDraw);
    m_vertexColor.bind();
    m_vertexColor.allocate(colors, sizeClouds* 3 * sizeof(GLfloat));
    m_vertexColor.release();

    if(!onlyColor)
        delete [] vertices;
    delete [] colors;
}

void GLWidget::setData(cData *data)
{
    m_Data = data;

    if (m_Data->NbClouds())
    {
        m_bDisplayMode2D = false;

        setBufferGl();
 
        setZoom(m_Data->getCloud(0)->getScale());

        m_params.m_translationMatrix[0] = -m_Data->m_cX;
        m_params.m_translationMatrix[1] = -m_Data->m_cY;
        m_params.m_translationMatrix[2] = -m_Data->m_cZ;
    }

    if (m_Data->NbImages())
    {
        m_bDisplayMode2D = true;

        _glImg = QGLWidget::convertToGLFormat( *m_Data->getCurImage() );

        applyGamma(m_params.getGamma());

        zoomFit();

        //position de l'image dans la vue gl             
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        glGetDoublev (GL_MODELVIEW_MATRIX, _mvmatrix);
        glGenTextures(1, &m_textureImage );

        ImageToTexture(m_textureImage, &_glImg);

        if(_mask)
            delete _mask;
        else
           glGenTextures(1, &m_textureMask );

        _mask = new QImage(_glImg.size(),_glImg.format());

        if (m_Data->NbMasks())
        {
            *_mask = QGLWidget::convertToGLFormat( *m_Data->getCurMask() );
            m_bFirstAction = false;
        }
        else
        {
            QGLWidget::convertToGLFormat(*_mask);
            _mask->fill(Qt::white);
            m_bFirstAction = true;
        }

        ImageToTexture(m_textureMask, _mask);
    }

    if (m_Data->NbCameras())
    {
        m_bDisplayMode2D = false;
        //TODO
    }

    glGetIntegerv (GL_VIEWPORT, _glViewport);
}

void GLWidget::dragEnterEvent(QDragEnterEvent *event)
{
    const QMimeData* mimeData = event->mimeData();

    if (mimeData->hasFormat("text/uri-list"))
        event->acceptProposedAction();
}

void GLWidget::dropEvent(QDropEvent *event)
{
    const QMimeData* mimeData = event->mimeData();

    if (mimeData->hasFormat("text/uri-list"))
    {
        QByteArray data = mimeData->data("text/uri-list");
        QStringList fileNames = QUrl::fromPercentEncoding(data).split(QRegExp("\\n+"),QString::SkipEmptyParts);

        for (int i=0;i<fileNames.size();++i)
        {
            fileNames[i] = fileNames[i].trimmed();

#if defined(_WIN32) || defined(WIN32)
            fileNames[i].remove("file:///");
#else
            fileNames[i].remove("file://");
#endif

#ifdef _DEBUG
            QString formatedMessage = QString("File dropped: %1").arg(fileNames[i]);
            printf(" %s\n",qPrintable(formatedMessage));
#endif
        }

        if (!fileNames.empty())
            emit filesDropped(fileNames);

        setFocus();

        event->acceptProposedAction();
    }

    event->ignore();
}

void GLWidget::displayNewMessage(const QString& message,
                                 MessagePosition pos)
{
    if (message.isEmpty())
    {
        m_messagesToDisplay.clear();

        return;
    }

    MessageToDisplay mess;
    mess.message = message;
    mess.position = pos;
    m_messagesToDisplay.push_back(mess);
}

void GLWidget::drawGradientBackground()
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    int w = (_glViewport[2]>>1)+1;
    int h = (_glViewport[3]>>1)+1;
    glOrtho(-w,w,-h,h,-2.f, 2.f);

    const uchar BkgColor[3] = {(uchar) colorBG0.red(),(uchar) colorBG0.green(), (uchar) colorBG0.blue()};
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    //Gradient "texture" drawing
    glBegin(GL_QUADS);
    //user-defined background color for gradient start
    glColor3ubv(BkgColor);
    glVertex2f(-w,h);
    glVertex2f(w,h);
    //and the inverse of points color for gradient end
    glColor3ub(colorBG1.red(),colorBG1.green(),colorBG1.blue());
    glVertex2f(w,-h);
    glVertex2f(-w,-h);
    glEnd();
}

void GLWidget::drawPolygon()
{
    glColor3f(.1f,1.f,.2f);

    glLineWidth(1.0f);

    if (m_Data->NbImages())
    {
        QVector <QPointF> poly;
        for (int aK = 0;aK < (int) m_polygon.size(); ++aK)
        {
            poly.push_back(ImageToWindow(m_polygon[aK]));
        }

        glBegin(m_bPolyIsClosed ? GL_LINE_LOOP : GL_LINE_STRIP);
        for (int aK = 0;aK < (int) poly.size(); ++aK)
        {
            glVertex2f(poly[aK].x(), poly[aK].y());
        }
        glEnd();

        glColor3f(1.f,0.f,0.f);

        if (_idx >=0)
        {
            for (int aK = 0;aK < _idx; ++aK)
                glDrawUnitCircle(2, poly[aK].x(), poly[aK].y(), 3.0, 8);

            glColor3f(0.f,0.f,1.f);
            glDrawUnitCircle(2, poly[_idx].x(), poly[_idx].y(), 3.0, 8);

            glColor3f(1.f,0.f,0.f);
            for (int aK = _idx+1;aK < (int) poly.size(); ++aK)
                glDrawUnitCircle(2, poly[aK].x(), poly[aK].y(), 3.0, 8);
        }
        else
        {
            for (int aK = 0;aK < (int) poly.size(); ++aK)
               glDrawUnitCircle(2, poly[aK].x(), poly[aK].y(), 3.0, 8);
        }
    }
    else
    {
        glBegin(m_bPolyIsClosed ? GL_LINE_LOOP : GL_LINE_STRIP);
        for (int aK = 0;aK < (int) m_polygon.size(); ++aK)
        {
            glVertex2f(m_polygon[aK].x(), m_polygon[aK].y());
        }
        glEnd();

        glColor3f(1.f,0.f,0.f);

        if (_idx >=0)
        {
            for (int aK = 0;aK < _idx; ++aK)
                glDrawUnitCircle(2, m_polygon[aK].x(), m_polygon[aK].y(), 3.0, 8);

            glColor3f(0.f,0.f,1.f);
            glDrawUnitCircle(2, m_polygon[_idx].x(), m_polygon[_idx].y(), 3.0, 8);

            glColor3f(1.f,0.f,0.f);
            for (int aK = _idx+1;aK < (int) m_polygon.size(); ++aK)
                glDrawUnitCircle(2, m_polygon[aK].x(), m_polygon[aK].y(), 3.0, 8);
        }
        else
        {
            for (int aK = 0;aK < (int) m_polygon.size(); ++aK)
               glDrawUnitCircle(2, m_polygon[aK].x(), m_polygon[aK].y(), 3.0, 8);
        }
    }
}

void GLWidget::drawPointAndSegments()
{
    if (m_polygon2.size())
    {
        glLineStipple(2, 0xAAAA);
        glEnable(GL_LINE_STIPPLE);

        glEnable (GL_BLEND);
        glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glHint (GL_LINE_SMOOTH_HINT, GL_DONT_CARE);

        glColor3f(0.f,0.f,1.f);

        glLineWidth(1.0f);

        glBegin(GL_LINE_STRIP);
        for (int aK=0;aK < (int) m_polygon2.size(); ++aK)
        {
            glVertex2f(m_polygon2[aK].x(), m_polygon2[aK].y());
        }
        glEnd();

        glDisable(GL_LINE_STIPPLE);

        enableOptionLine();

        glDrawUnitCircle(2, m_polygon2[1].x(), m_polygon2[1].y(), 3.0, 8);

        disableOptionLine();
    }
}

// zoom in 3D mode
void GLWidget::zoom()
{
    GLdouble zoom = m_params.zoom;

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    glOrtho(-zoom*m_glRatio,zoom*m_glRatio,-zoom, zoom,-2.f*m_Data->m_diam, 2.f*m_Data->m_diam);

    glMatrixMode(GL_MODELVIEW);
}

void GLWidget::setInteractionMode(INTERACTION_MODE mode)
{
    m_interactionMode = mode;

    switch (mode)
    {
    case TRANSFORM_CAMERA:
        setMouseTracking(false);
        removeEventFilter(this);
        break;
    case SELECTION:
        if(!m_Data->NbImages())
            setProjectionMatrix();
        installEventFilter(this);
        setMouseTracking(true);
        break;
    default:
        break;
    }
}

void GLWidget::setView(VIEW_ORIENTATION orientation)
{
    makeCurrent();

    GLdouble eye[3] = {0.0, 0.0, 0.0};
    GLdouble top[3] = {0.0, 0.0, 0.0};
    GLdouble s[3]   = {0.0, 0.0, 0.0};
    GLdouble u[3]   = {0.0, 0.0, 0.0};

    switch (orientation)
    {
    case TOP_VIEW:
        eye[2] = -1.0;
        top[1] =  1.0;
        break;
    case BOTTOM_VIEW:
        eye[2] =  1.0;
        top[1] = -1.0;
        break;
    case FRONT_VIEW:
        eye[1] = 1.0;
        top[2] = 1.0;
        break;
    case BACK_VIEW:
        eye[1] = -1.0;
        top[2] =  1.0;
        break;
    case LEFT_VIEW:
        eye[0] = 1.0;
        top[2] = 1.0;
        break;
    case RIGHT_VIEW:
        eye[0] = -1.0;
        top[2] =  1.0;
    }

    crossprod(eye, top, s);
    crossprod(s, eye, u);

    _m_g_rotationMatrix[0] = s[0];
    _m_g_rotationMatrix[1] = s[1];
    _m_g_rotationMatrix[2] = s[2];

    _m_g_rotationMatrix[3] = u[0];
    _m_g_rotationMatrix[4] = u[1];
    _m_g_rotationMatrix[5] = u[2];

    _m_g_rotationMatrix[6] = -eye[0];
    _m_g_rotationMatrix[7] = -eye[1];
    _m_g_rotationMatrix[8] = -eye[2];

    m_params.m_translationMatrix[0] = m_Data->m_cX;
    m_params.m_translationMatrix[1] = m_Data->m_cY;
    m_params.m_translationMatrix[2] = m_Data->m_cZ;
}

void GLWidget::onWheelEvent(float wheelDelta_deg)
{
    //convert degrees in zoom 'power'
    float zoomFactor = pow(1.1f,wheelDelta_deg *.05f);

    setZoom(m_params.zoom*zoomFactor);
}

void GLWidget::setZoom(float value)
{
    if (value < GL_MIN_ZOOM_RATIO)
        value = GL_MIN_ZOOM_RATIO;
    else if (value > GL_MAX_ZOOM_RATIO)
        value = GL_MAX_ZOOM_RATIO;

    m_params.zoom = value;   

    update();
}

void GLWidget::zoomFit()
{
    //width and height ratio between viewport and image
    m_rw = (float)_glImg.width()/ _glViewport[2];
    m_rh = (float)_glImg.height()/_glViewport[3];

    if(m_rw>m_rh)
        setZoom(1.f/m_rw); //orientation landscape
    else
        setZoom(1.f/m_rh); //orientation portrait

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glPushMatrix();
    glScalef(m_params.zoom, m_params.zoom, 1.0);
    glTranslatef(-m_rw,-m_rh,0);
    glGetDoublev (GL_PROJECTION_MATRIX, _projmatrix);
    glPopMatrix();

    m_glPosition[0] = 0;
    m_glPosition[1] = 0;

}

void GLWidget::zoomFactor(int percent)
{
    setZoom((float) percent / 100.f);
}

void GLWidget::wheelEvent(QWheelEvent* event)
{
    if (m_interactionMode == SELECTION)
    {
        event->ignore();
        return;
    }

    //see QWheelEvent documentation ("distance that the wheel is rotated, in eighths of a degree")
    float wheelDelta_deg = (float)event->delta() / 8.f;

    _m_lastPosZoom = event->pos();

    onWheelEvent(wheelDelta_deg);
}

void GLWidget::ImageToTexture(GLuint idTexture, QImage *image)
{
    glBindTexture( GL_TEXTURE_2D, idTexture );
    glTexImage2D( GL_TEXTURE_2D, 0, 4, image->width(), image->height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, image->bits());
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
    glBindTexture( GL_TEXTURE_2D, 0);
    glDisable(GL_TEXTURE_2D);
}

void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
    if (event->x()<0 || event->y()<0 || event->x()>width() || event->y()>height())
        return;

    QPoint pos = event->pos();

    if (m_interactionMode == SELECTION)
    {
        findClosestPoint();
        update();
    }
    else
    {
        QPoint dp = pos-m_lastPos;

        if ( _m_g_mouseLeftDown ) // rotation autour de X et Y
        {
            float d_angleX = m_speed * dp.y() / (float) _glViewport[3];
            float d_angleY = m_speed * dp.x() / (float) _glViewport[2];

            m_params.angleX += d_angleX;
            m_params.angleY += d_angleY;

            setRotateOx_m33( d_angleX, _m_g_rotationOx );
            setRotateOy_m33( d_angleY, _m_g_rotationOy );

            mult_m33( _m_g_rotationOx, _m_g_rotationMatrix, _m_g_tmpoMatrix );
            mult_m33( _m_g_rotationOy, _m_g_tmpoMatrix, _m_g_rotationMatrix );
        }
        else if ( _m_g_mouseMiddleDown )
        {
            if (event->modifiers() & Qt::ShiftModifier) // zoom
            {
                if (dp.y() > 0) m_params.zoom *= pow(2.f, dp.y() *.05f);
                else if (dp.y() < 0) m_params.zoom /= pow(2.f, -dp.y() *.05f);
            }
            else if((_glViewport[2]!=0) || (_glViewport[3]!=0)) // translation
            {
                    if (m_Data->NbImages())
                    {
                        m_glPosition[0] += 2.0f*( (float)dp.x()/(_glViewport[2]*m_params.zoom) );
                        m_glPosition[1] -= 2.0f*( (float)dp.y()/(_glViewport[3]*m_params.zoom) );
                    }
                    else
                    {
                        m_bObjectCenteredView = false;
                        m_params.m_translationMatrix[0] += m_speed * dp.x()*m_Data->m_diam/_glViewport[2];
                        m_params.m_translationMatrix[1] -= m_speed * dp.y()*m_Data->m_diam/_glViewport[3];
                    }
            }
        }
        else if ( _m_g_mouseRightDown ) // rotation autour de Z
        {
            float d_angleZ =  m_speed * dp.x() / (float) _glViewport[2];

            m_params.angleZ += d_angleZ;

            setRotateOz_m33( d_angleZ, _m_g_rotationOz );

            mult_m33( _m_g_rotationOz, _m_g_rotationMatrix, _m_g_tmpoMatrix );

            for (int i = 0; i < 9; ++i) _m_g_rotationMatrix[i] = _m_g_tmpoMatrix[i];
        }

        update();
    }

    m_lastPos = pos;
}

void GLWidget::mouseDoubleClickEvent(QMouseEvent *event)
{
    if (m_Data->NbClouds())
    {
        QPointF pos = event->localPos();

        setProjectionMatrix();

        int idx1 = -1;
        int idx2;

        pos.setY(_glViewport[3] - pos.y());

        for (int aK=0; aK < m_Data->NbClouds();++aK)
        {
            float sqrD;
            float dist = FLT_MAX;
            idx2 = -1;
            QPointF proj;

            Cloud *a_cloud = m_Data->getCloud(aK);

            for (int bK=0; bK < a_cloud->size();++bK)
            {
                getProjection(proj, a_cloud->getVertex( bK ));

                sqrD = (proj.x()-pos.x())*(proj.x()-pos.x()) + (proj.y()-pos.y())*(proj.y()-pos.y());

                if (sqrD < dist )
                {
                    dist = sqrD;
                    idx1 = aK;
                    idx2 = bK;
                }
            }
        }

        if ((idx1>=0) && (idx2>=0))
        {
            //final center:
            Cloud *a_cloud = m_Data->getCloud(idx1);
            Vertex &P = a_cloud->getVertex( idx2 );

            m_Data->m_cX = P.x();
            m_Data->m_cY = P.y();
            m_Data->m_cZ = P.z();

            update();
        }
    }
}

bool isPointInsidePoly(const QPointF& P, const QVector< QPointF> poly)
{
    int vertices=poly.size();
    if (vertices<3)
        return false;

    bool inside = false;

    QPointF A = poly[0];
    for (int i=1;i<=vertices;++i)
    {
        QPointF B = poly[i%vertices];

        //Point Inclusion in Polygon Test (inspired from W. Randolph Franklin - WRF)
        if (((B.y() <= P.y()) && (P.y()<A.y())) ||
                ((A.y() <= P.y()) && (P.y()<B.y())))
        {
            float ABy = A.y()-B.y();
            float t = (P.x()-B.x())*ABy-(A.x()-B.x())*(P.y()-B.y());
            if (ABy<0)
                t=-t;

            if (t<0)
                inside = !inside;
        }

        A=B;
    }

    return inside;
}

void GLWidget::setProjectionMatrix()
{
    glMatrixMode(GL_MODELVIEW);
    glGetDoublev(GL_MODELVIEW_MATRIX, _mvmatrix);

    glMatrixMode(GL_PROJECTION);
    glGetDoublev(GL_PROJECTION_MATRIX, _projmatrix);

    glGetIntegerv(GL_VIEWPORT, _glViewport);
}

void GLWidget::getProjection(QPointF &P2D, Vertex P)
{
    GLdouble xp,yp,zp;
    gluProject(P.x(),P.y(),P.z(),_mvmatrix,_projmatrix,_glViewport,&xp,&yp,&zp);
    P2D = QPointF(xp,yp);
}

QPointF GLWidget::WindowToImage(QPointF const &pt)
{
    QPointF res;

    res.setX(( pt.x()         - _glViewport[2]*.5f ) - _projmatrix[12]*_glViewport[2]*.5f);
    res.setY((-pt.y()  -1.f   + _glViewport[3]*.5f ) - _projmatrix[13]*_glViewport[3]*.5f);

    res /= m_params.zoom;

    return res;
}

QPointF GLWidget::ImageToWindow(QPointF const &im)
{
    QPointF res;

    res.setX(im.x()*m_params.zoom + _glViewport[2]*.5f + _projmatrix[12]*_glViewport[2]*.5f);
    res.setY(- 1.f - im.y()*m_params.zoom + _glViewport[3]*.5f - _projmatrix[13]*_glViewport[3]*.5f);

    return res;
}

void GLWidget::Select(int mode)
{
    QPointF P2D;
    bool pointInside;
    QVector< QPointF> polyg;

    if(mode == ADD || mode == SUB)
    {
        if ((m_polygon.size() < 3) || (!m_bPolyIsClosed))
            return;

        if (!m_bDisplayMode2D)
        {
            for (int aK=0; aK < (int) m_polygon.size(); ++aK)
            {
               polyg.push_back(QPointF(m_polygon[aK].x(), _glViewport[3] - m_polygon[aK].y()));
            }
        }
        else
        {
            for (int aK=0; aK < (int) m_polygon.size(); ++aK)
            {
                polyg.push_back(m_polygon[aK]);
            }
        }

    }

    if (m_bDisplayMode2D)
    {
         QPainter    p;
         QBrush SBrush(Qt::white);
         QBrush NSBrush(Qt::black);

         p.begin(_mask);
         p.setCompositionMode(QPainter::CompositionMode_Source);
         p.setPen(Qt::NoPen);

         if(mode == ADD)
         {
             if (m_bFirstAction)
             {
                 p.fillRect(_mask->rect(), Qt::black);
             }
             p.setBrush(SBrush);
             p.drawPolygon(polyg.data(),polyg.size());
         }
         else if(mode == SUB)
         {
             p.setBrush(NSBrush);
             p.drawPolygon(polyg.data(),polyg.size());
         }
         else if(mode == ALL)
         {
             p.fillRect(_mask->rect(), Qt::white);
         }
         else if(mode == NONE)
         {
             p.fillRect(_mask->rect(), Qt::black);
         }
         p.end();

         if(mode == INVERT)         
            _mask->invertPixels(QImage::InvertRgb);

         ImageToTexture(m_textureMask, _mask);
    }
    else
    {
        for (int aK=0; aK < m_Data->NbClouds(); ++aK)
        {
            Cloud *a_cloud = m_Data->getCloud(aK);

            for (int bK=0; bK < a_cloud->size();++bK)
            {
                Vertex &P = a_cloud->getVertex( bK );
                switch (mode)
                {
                case ADD:
                    getProjection(P2D, P);
                    pointInside = isPointInsidePoly(P2D,polyg);
                    if (m_bFirstAction)
                        emit selectedPoint((uint)aK,(uint)bK,pointInside);
                    else
                        emit selectedPoint((uint)aK,(uint)bK,pointInside||P.isVisible());
                    break;
                case SUB:
                    if (P.isVisible())
                    {
                        getProjection(P2D, P);
                        pointInside = isPointInsidePoly(P2D,polyg);
                        emit selectedPoint((uint)aK,(uint)bK,!pointInside);
                    }
                    break;
                case INVERT:
                    //if (m_previousAction == NONE)  m_bFirstAction = true;
                    emit selectedPoint((uint)aK,(uint)bK,!P.isVisible());
                    break;
                case ALL:
                    m_bFirstAction = true;
                    emit selectedPoint((uint)aK,(uint)bK, true);
                    break;
                case NONE:
                    emit selectedPoint((uint)aK,(uint)bK,false);
                    break;
                }
            }

            setBufferGl(true);
        }
    }

    if (((mode == ADD)||(mode == SUB)) && (m_bFirstAction)) m_bFirstAction = false;

    selectInfos info;
    info.params = m_params;
    info.poly   = m_polygon;
    info.selection_mode   = mode;

    m_infos.push_back(info);

    clearPolyline();
}

void GLWidget::clearPolyline()
{
    m_polygon.clear();
    m_polygon2.clear();
    m_Click = 0;
    _idx = -1;
    m_bPolyIsClosed = false;
    update();
}

void GLWidget::closePolyline()
{
    if (!m_bPolyIsClosed)
    {
        //remove last point if needed
        int sz = m_polygon.size();
        if (sz > 2) m_polygon.resize(sz-1);

        m_bPolyIsClosed = true;
    }
}

void GLWidget::ptSizeUp(bool up)
{
    if (up)
        m_params.PointSize++;
    else
        m_params.PointSize--;

    if (m_params.PointSize == 0)
        m_params.PointSize = 1;

    glPointSize(m_params.PointSize);

    update();
}

void GLWidget::drawAxis()
{
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();

    if (m_trihedronGLList == GL_INVALID_LIST_ID)
    {
        m_trihedronGLList = glGenLists(1);
        glNewList(m_trihedronGLList, GL_COMPILE);

        glPushAttrib(GL_LINE_BIT | GL_DEPTH_BUFFER_BIT);
        glEnable(GL_LINE_SMOOTH);
        glEnable(GL_DEPTH_TEST);

        //trihedra OpenGL drawing
        glBegin(GL_LINES);
        glColor3f(1.0f,0.0f,0.0f);
        glVertex3f(0.0f,0.0f,0.0f);
        glVertex3f(0.4f,0.0f,0.0f);
        glColor3f(0.0f,1.0f,0.0f);
        glVertex3f(0.0f,0.0f,0.0f);
        glVertex3f(0.0f,0.4f,0.0f);
        glColor3f(0.0f,0.7f,1.0f);
        glVertex3f(0.0f,0.0f,0.0f);
        glVertex3f(0.0f,0.0f,0.4f);
        glEnd();

        glPopAttrib();

        glEndList();
    }

    glTranslatef(m_Data->m_cX,m_Data->m_cY,m_Data->m_cZ);

    glCallList(m_trihedronGLList);

    glPopMatrix();
}

void GLWidget::drawBall()
{
    if (!m_bObjectCenteredView) return;

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();    

    // ball radius
    float scale = m_Data->m_diam / 1.5f;

    if (m_ballGLList == GL_INVALID_LIST_ID)
    {
        incrNbGLLists();
        m_ballGLList = getNbGLLists();
        glNewList(m_ballGLList, GL_COMPILE);

        //draw 3 circles
        glPushAttrib(GL_LINE_BIT);
        glEnable(GL_LINE_SMOOTH);
        glPushAttrib(GL_COLOR_BUFFER_BIT);
        glEnable(GL_BLEND);
        const float c_alpha = 0.6f;
        glLineWidth(1.0f);

        glColor4f(1.0f,0.0f,0.0f,c_alpha);
        glDrawUnitCircle(0, 0, 0, 1.f);
        glBegin(GL_LINES);
        glVertex3f(-1.0f,0.0f,0.0f);
        glVertex3f( 1.0f,0.0f,0.0f);
        glEnd();

        glColor4f(0.0f,1.0f,0.0f,c_alpha);
        glDrawUnitCircle(1, 0, 0, 1.f);
        glBegin(GL_LINES);
        glVertex3f(0.0f,-1.0f,0.0f);
        glVertex3f(0.0f, 1.0f,0.0f);
        glEnd();

        glColor4f(0.0f,0.7f,1.0f,c_alpha);
        glDrawUnitCircle(2, 0, 0, 1.f);
        glBegin(GL_LINES);
        glVertex3f(0.0f,0.0f,-1.0f);
        glVertex3f(0.0f,0.0f, 1.0f);
        glEnd();

        glPopAttrib();

        glEndList();
    }

    glTranslatef(m_Data->m_cX,m_Data->m_cY,m_Data->m_cZ);
    glScalef(scale,scale,scale);

    glCallList(m_ballGLList);

    glPopMatrix();
    glDisable(GL_BLEND);
}

void GLWidget::drawCams()
{
    float scale = .1f;

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();

    incrNbGLLists();
    GLuint list = getNbGLLists();
    glNewList(list, GL_COMPILE);

    glLineWidth(2);
    glPointSize(7);
    if (m_Data->NbClouds())
        scale = .1f*m_Data->getCloud(0)->getScale();

    for (int i=0; i<m_Data->NbCameras();i++)
    {
        CamStenope * pCam = m_Data->getCamera(i);

        //REAL f = pCam->Focale();
        Pt3dr C  = pCam->VraiOpticalCenter();
        Pt3dr P1 = pCam->ImEtProf2Terrain(Pt2dr(0,0),scale);
        Pt3dr P2 = pCam->ImEtProf2Terrain(Pt2dr(pCam->Sz().x,0),scale);
        Pt3dr P3 = pCam->ImEtProf2Terrain(Pt2dr(0,pCam->Sz().y),scale);
        Pt3dr P4 = pCam->ImEtProf2Terrain(Pt2dr(pCam->Sz().x,pCam->Sz().y),scale);

        //translation
        if (m_Data->NbClouds())
        {
            Pt3dr translation = m_Data->getCloud(0)->getTranslation();

            C = C + translation;
            P1 = P1 + translation;
            P2 = P2 + translation;
            P3 = P3 + translation;
            P4 = P4 + translation;
        }

        glBegin(GL_LINES);
            //perspective cone
            qglColor(QColor(0,0,0));
            glVertex3d(C.x, C.y, C.z);
            glVertex3d(P1.x, P1.y, P1.z);

            glVertex3d(C.x, C.y, C.z);
            glVertex3d(P2.x, P2.y, P2.z);

            glVertex3d(C.x, C.y, C.z);
            glVertex3d(P3.x, P3.y, P3.z);

            glVertex3d(C.x, C.y, C.z);
            glVertex3d(P4.x, P4.y, P4.z);

            //Image
            qglColor(QColor(255,0,0));
            glVertex3d(P1.x, P1.y, P1.z);
            glVertex3d(P2.x, P2.y, P2.z);

            glVertex3d(P4.x, P4.y, P4.z);
            glVertex3d(P2.x, P2.y, P2.z);

            glVertex3d(P3.x, P3.y, P3.z);
            glVertex3d(P1.x, P1.y, P1.z);

            glVertex3d(P4.x, P4.y, P4.z);
            glVertex3d(P3.x, P3.y, P3.z);
        glEnd();

        glBegin(GL_POINTS);
            glVertex3d(C.x, C.y, C.z);
        glEnd();
    }

    glEndList();

    glCallList(list);

    glPointSize(m_params.PointSize);
    glLineWidth(m_params.LineWidth);
    glPopMatrix();
}

void GLWidget::drawBbox()
{
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();

    incrNbGLLists();
    GLuint list = getNbGLLists();
    glNewList(list, GL_COMPILE);

    glLineWidth(1);

    if (m_Data->NbClouds())
    {
        float minX, minY, minZ, maxX, maxY, maxZ;

        minX = m_Data->m_minX;
        minY = m_Data->m_minY;
        minZ = m_Data->m_minZ;
        maxX = m_Data->m_maxX;
        maxY = m_Data->m_maxY;
        maxZ = m_Data->m_maxZ;

        Pt3dr P1(minX, minY, minZ);
        Pt3dr P2(minX, minY, maxZ);
        Pt3dr P3(minX, maxY, maxZ);
        Pt3dr P4(minX, maxY, minZ);
        Pt3dr P5(maxX, minY, minZ);
        Pt3dr P6(maxX, maxY, minZ);
        Pt3dr P7(maxX, maxY, maxZ);
        Pt3dr P8(maxX, minY, maxZ);

        glBegin(GL_LINES);

        qglColor(QColor("orange"));

            glVertex3d(P1.x, P1.y, P1.z);
            glVertex3d(P2.x, P2.y, P2.z);

            glVertex3d(P3.x, P3.y, P3.z);
            glVertex3d(P2.x, P2.y, P2.z);

            glVertex3d(P1.x, P1.y, P1.z);
            glVertex3d(P4.x, P4.y, P4.z);

            glVertex3d(P1.x, P1.y, P1.z);
            glVertex3d(P5.x, P5.y, P5.z);

            glVertex3d(P7.x, P7.y, P7.z);
            glVertex3d(P3.x, P3.y, P3.z);

            glVertex3d(P7.x, P7.y, P7.z);
            glVertex3d(P6.x, P6.y, P6.z);

            glVertex3d(P8.x, P8.y, P8.z);
            glVertex3d(P5.x, P5.y, P5.z);

            glVertex3d(P7.x, P7.y, P7.z);
            glVertex3d(P8.x, P8.y, P8.z);

            glVertex3d(P5.x, P5.y, P5.z);
            glVertex3d(P6.x, P6.y, P6.z);

            glVertex3d(P4.x, P4.y, P4.z);
            glVertex3d(P6.x, P6.y, P6.z);

            glVertex3d(P8.x, P8.y, P8.z);
            glVertex3d(P2.x, P2.y, P2.z);

            glVertex3d(P4.x, P4.y, P4.z);
            glVertex3d(P3.x, P3.y, P3.z);

        glEnd();

        glEndList();
    }

    glCallList(list);

    glPointSize(m_params.PointSize);
    glLineWidth(m_params.LineWidth);
    glPopMatrix();
}

void GLWidget::showAxis(bool show)
{
    m_bDrawAxis = show;
    if (m_bDrawAxis) m_bDrawBall = false;

    update();
}

void GLWidget::showBall(bool show)
{
    m_bDrawBall = show;
    if (m_bDrawBall) m_bDrawAxis = false;

    update();
}

void GLWidget::showCams(bool show)
{
    m_bDrawCams = show;

    update();
}

void GLWidget::showBBox(bool show)
{
    m_bDrawBbox = show;

    update();
}

void GLWidget::showMessages(bool show)
{
    m_bDrawMessages = show;

    if (m_bDrawMessages)
    {
        if (m_interactionMode == TRANSFORM_CAMERA)
            showMoveMessages();
        else if (m_interactionMode == SELECTION)
            showSelectionMessages();
    }
    else
        displayNewMessage(QString());

    update();
}

bool GLWidget::showMessages(){return m_bDrawMessages;}

void GLWidget::showSelectionMessages()
{
    displayNewMessage(QString());
    displayNewMessage(tr("Selection mode"),UPPER_CENTER_MESSAGE);
    displayNewMessage(tr("Left click: add contour point / Right click: close"),LOWER_CENTER_MESSAGE);
    displayNewMessage(tr("Space: add / Suppr: delete"),LOWER_CENTER_MESSAGE);
}

void GLWidget::showMoveMessages()
{
    displayNewMessage(QString());
    displayNewMessage(tr("Move mode"),UPPER_CENTER_MESSAGE);
    if (m_bDisplayMode2D)
        displayNewMessage(tr("Wheel: zoom / Right click: translate viewpoint"),LOWER_CENTER_MESSAGE);
    else
        displayNewMessage(tr("Left click: rotate viewpoint / Right click: translate viewpoint"),LOWER_CENTER_MESSAGE);
}

void GLWidget::reset()
{
    _m_g_rotationMatrix[0] = _m_g_rotationMatrix[4] = _m_g_rotationMatrix[8] = 1;
    _m_g_rotationMatrix[1] = _m_g_rotationMatrix[2] = _m_g_rotationMatrix[3] = 0;
    _m_g_rotationMatrix[5] = _m_g_rotationMatrix[6] = _m_g_rotationMatrix[7] = 0;

    clearPolyline();

    m_params.reset();
    m_Data->clearClouds();
    m_Data->clearCameras();
    m_Data->clearImages();
    m_Data->clearMasks();

    m_bFirstAction = true;
}

void GLWidget::applyGamma(float aGamma)
{
    if (aGamma == 1.f) return;

    QRgb  pixel;
    int r,g,b;

    for(int i=0; i< _glImg.width();++i)
        for(int j=0; j<_glImg.height();++j)
        {
            pixel = _glImg.pixel(i,j);

            r = 255*pow((float) qRed(pixel)  / 255.f, 1.f / aGamma);
            g = 255*pow((float) qGreen(pixel)/ 255.f, 1.f / aGamma);
            b = 255*pow((float) qBlue(pixel) / 255.f, 1.f / aGamma);

            if (r>255) r = 255;
            if (g>255) g = 255;
            if (b>255) b = 255;

            _glImg.setPixel(i,j, qRgb(r,g,b) );
        }
}

float segmentDistToPoint(QPointF segA, QPointF segB, QPointF p)
{
    QPointF p2(segB.x() - segA.x(), segB.y() - segA.y());
    float nrm = (p2.x()*p2.x() + p2.y()*p2.y());
    float u = ((p.x() - segA.x()) * p2.x() + (p.y() - segA.y()) * p2.y()) / nrm;

    if (u > 1)
        u = 1;
    else if (u < 0)
        u = 0;

    float x = segA.x() + u * p2.x();
    float y = segA.y() + u * p2.y();

    float dx = x - p.x();
    float dy = y - p.y();

    return sqrt(dx*dx + dy*dy);
}

void GLWidget::fillPolygon2()
{
    float dist, dist2;
    dist2 = FLT_MAX;
    int idx = -1;

    QVector < QPointF > polygon = m_polygon;
    polygon.push_back(polygon[0]);

    for (int aK =0; aK < (int) polygon.size()-1;++aK)
    {
        dist = segmentDistToPoint(polygon[aK], polygon[aK+1], m_lastPos);

        if (dist < dist2)
        {
            dist2 = dist;
            idx = aK;
        }
    }

    if (idx != -1)
    {
        m_polygon2.clear();
        m_polygon2.push_back(polygon[idx]);
        m_polygon2.push_back(m_lastPos);
        m_polygon2.push_back(polygon[idx+1]);
    }
}

void GLWidget::findClosestPoint()
{
    _idx = -1;
    float dist, dist2;
    dist2 = (float) m_sqr_radius;

    for (int aK = 0; aK < (int) m_polygon.size();++aK)
    {
        dist  = (float)(m_lastPos.x() - m_polygon[aK].x())*(m_lastPos.x() - m_polygon[aK].x()) +
                (m_lastPos.y() - m_polygon[aK].y())*(m_lastPos.y() - m_polygon[aK].y());

        if  (dist < dist2)
        {
            dist2 = dist;
            _idx = aK;
        }
    } 
}
