#include <fstream>
#include <string>
#include <iostream>
#include <math.h>
#include <memory>
#include <iomanip>
#include <vector>
#include <unistd.h>


typedef enum class enm
{
    FAILEDEXITCOMMAND = 1, COMMANDNOTFOUND, SETPATTERN
}ERRORCODE;

template<int n>
class Vec
{
    public:
    int static const dimension = n;
    double a[n];
    constexpr Vec()
    {
        for(int i=0; i<dimension; ++i)
            a[i] = 0;
    }
    constexpr Vec(double* r)
    {
        for(int i=0; i<dimension; ++i)
            a[i] = r[i];
    }
    constexpr double length() const
    {   
        double sum = 0;
        for(int i=0; i<dimension; ++i)
            sum += a[i]*a[i];
        return sqrt(sum);
    }
    constexpr Vec normalize();
};

template<int n>
std::ostream& operator<<(std::ostream& l, const Vec<n>& r) 
{
    std::cout << '(';
    for(int i=0; i<r.dimension-1; i++)
        std::cout << r.a[i] << ", ";
    std::cout << r.a[r.dimension-1];
    std::cout << ')';
    return l;
}

template<int n>
constexpr Vec<n> operator / (Vec<n> const& l, double const& r)
{
    Vec<n> ret;
    for(int i=0; i<l.dimension; i++)
        ret.a[i] = l.a[i] / r;
    return ret;
}

template<int n>
constexpr void operator /=(Vec<n>& l, double const& r)
{   l = l/r;}

template<int n>
constexpr Vec<n> Vec<n>::normalize()
{   *this /= this->length();    return *this;}

template<int n>
constexpr Vec<n> operator - (Vec<n> const& l, Vec<n> const& r)
{   
    Vec<n> ret;
    for(int i=0; i<l.dimension; i++)
        ret.a[i] = l.a[i] - r.a[i];
    return ret;
}

template<int n>
constexpr Vec<n> operator * (Vec<n> const& l, double const& r)
{
    Vec<n> ret;
    for(int i=0; i<l.dimension; i++)
        ret.a[i] = l.a[i] * r;
    return ret;
}


template<int n>
constexpr void operator *=(Vec<n>& l, double const& r)
{   l = l*r;}

template<int n>
constexpr double dot(Vec<n> const& l, Vec<n> const& r)
{   
    double ret = 0;
    for(int i=0; i<l.dimension; i++)
        ret += l.a[i] * r.a[i];
    return ret;
}

template<int n>
std::istream& operator>>(std::istream& l, Vec<n>& r)
{
    for(int i=0; i<r.dimension; i++)
        std::cin >> r.a[i];
    return l;
}


class Vec2 : public Vec<2>
{
    public:
    double& x = a[0],  y = a[1]; 

    Vec2(){}
    constexpr Vec2(double const& x_, double const& y_)
    {
        x = x_;
        y = y_;
    }

    constexpr Vec2(Vec<2> const& v)
    {
        x = v.a[0];
        y = v.a[1];
    }
    
    constexpr Vec2& operator=(const Vec<2> & r)
    {
        x = r.a[0];   y = r.a[1];
        return *this;
    }

    constexpr Vec2& operator=(const Vec2 & r)
    {
        x = r.a[0];   y = r.a[1];
        return *this;
    }
};

class Complex
{
    public:
    double re;
    double im;
    constexpr Complex(): re(0), im(0){}
    constexpr Complex(double const& re_, double const& im_): re(re_), im(im_){}

    Complex operator ~()
    {
        return Complex(re, -im);
    }
};

constexpr double abs(Complex const& comp)
{
    return sqrt(comp.re * comp.re + comp.im * comp.im);
}

std::ostream& operator <<(std::ostream& l, Complex r)
{
    return l << r.re << std::showpos << r.im << 'i';
}

constexpr Complex operator +(Complex const& l, Complex const& r)
{
    return Complex(l.re + r.re,   l.im + r.im);
}

constexpr Complex operator *(Complex const& l, Complex const& r)
{
    return Complex(l.re * r.re  -  l.im * r.im,       l.re * r.im  +  l.im * r.re);
}

std::istream& operator>>(std::istream& l, Complex& r)
{
    return l >> r.re >> r.im;
}

struct RGB
{
    char r;
    char g;
    char b;
};

/*cinでRGB型に書き込めるように*/
std::istream& operator>>(std::istream& l, RGB& rgb)
{
    /*cinでchar型をそのまま使うと数値として受け取ってもらえない*/
    short r, g, b;
    std::cin >> r >> g >> b;
    rgb = {(char)r, (char)g, (char)b};
    return l;
}

/*coutでRGB型を表示できるように*/
std::ostream& operator<<(std::ostream& left, RGB const& right)
{
    return left << '(' << +(unsigned char)right.r << ", " << +(unsigned char)right.g << ", " << +(unsigned char)right.b << ')';
}

constexpr RGB operator*(RGB const& left, double const& right)
{
    return {char((unsigned char)left.r * right), char((unsigned char)left.g * right), char((unsigned char)left.b * right)};
}

/*convertコマンドがなかった場合に実行される*/
void cmdnotfound(std::string const& command, std::vector<std::string> const& howto, std::vector<std::string> const& environment)
{
    std::cerr << '\n' << "Command '" + command + "' not found, but can be installed with:" << "\n\n"; 
    for(size_t i=0; i<howto.size(); i++)
        std::cerr << howto[i]  << "\t\t\t\t#" << environment[i] << '\n';
    std::cerr << std::endl;
}

/*pngファイルに書き込む*/
/*filename: 保存するファイルの名前,  width: 画像の横幅,  height: 画像の縦幅,  color: 色配列の先頭ポインタ*/
void writePNG(std::string filename, int width, int height, RGB* color)
{
    /*コマンドの実行でエラーが発生したら*/
    #define CHECK_ERROR(a, x)                                       \
    if((a == -1) || (a == 127) || (WIFEXITED(a) == 0))              \
    {                                                               \
        std::cerr << "An error occured!" << std::endl;              \
        exit(x);                                                    \
    }

    int res;

    res = system("type convert >/dev/null 2>&1");
    CHECK_ERROR(res, (int)ERRORCODE::FAILEDEXITCOMMAND)
    if(WEXITSTATUS(res) != 0)
    {
        std::vector<std::string> howto = {"sudo apt install imagemagick", "brew install imagemagick"};
        std::vector<std::string> environment = {"Linux", "Mac"};
        cmdnotfound("convert", howto, environment);
        exit((int)ERRORCODE::COMMANDNOTFOUND);
    }

    std::ofstream file(filename+".ppm");
    file << "P6" << std::endl;
    file << width << " " << height << std::endl;
    file << "255" << std::endl;

    for (int i = 0; i < width*height; i++)
        file << color[i].r << color[i].g << color[i].b;
        
    file.close();

    res = system(("convert " + filename + ".ppm " + filename + ".png >/dev/null 2>&1").c_str());
    CHECK_ERROR(res, (int)ERRORCODE::FAILEDEXITCOMMAND)
    if(WEXITSTATUS(res) != 0)
        std::cout << "ppm→png変換失敗" << std::endl;
    

    res = system(("rm " + filename + ".ppm >/dev/null 2>&1").c_str());
    CHECK_ERROR(res, (int)ERRORCODE::FAILEDEXITCOMMAND)
    if(WEXITSTATUS(res) != 0)
        std::cout << "ppmファイルの削除失敗" << std::endl;
    
}



/*文字を表示しながら入力値を返す*/
template<typename T>
T input(std::string info)
{
    /*バッファクリアとエラーフラグクリアを行う*/
    #define CIN_RESET  std::cin.clear(); std::cin.ignore(256,'\n');

    /*文字列aをb秒間だけ表示する(改行しない)*/
    #define SHOW_ERROR(a, b) \
    std::cout <<"\r                                                \r"; \
    std::cout << a << std::flush;                                       \
    sleep(b);                                                           \
    std::cout <<"\r                                                \r";

    T variable;
    while(1)
    {
        std::cout << info << std::endl;
        //CIN_RESET; // バッファクリア + エラーフラグクリア
        std::cin >> variable;

        /*正常に入力されれば*/
        if(std::cin.good() == true)
        {  
            std::cout << variable << "が入力されました。よろしければ5秒以内にEnterを入力してください。";

            time_t t0, t;   /*時刻を格納する変数*/
            time(&t0);      /*初期時刻*/
            CIN_RESET;      /*バッファクリア + エラーフラグクリア*/
            std::cin.get(); /*何か入力されるのを待つ*/
            time(&t);       /*入力された時刻*/

            int static constexpr limit = 5; /*制限時間*/

            /*制限時間以内に入力されていれば*/
            if(t-t0 < limit)
                return variable;
            /*制限時間を超過していれば*/
            else 
            {
                SHOW_ERROR("タイムオーバーです。再入力してください。" , 2);
                continue;
            }
        }
        /*入力で異常が発生すれば*/
        else
        {
            CIN_RESET; // バッファクリア + エラーフラグクリア
            SHOW_ERROR("入力で異常発生。再入力してください。", 2);
            continue;
        }
    }
}



class Image
{
    public:
    int const width, height;

    //width, heightの値でメモリの量が決まるのでdataを下で宣言すること

    private:
    RGB** data;

    public:

    /*コピーとムーブの禁止*/
    Image(Image const&) = delete;
    Image(Image const&&) = delete;
    Image& operator=(Image const&) = delete;
    Image& operator=(Image const&&) = delete;


    Image(int width_, int height_): width(width_), height(height_)
    {
        /*2次元配列をnewで確保したい*/
        /*そのデータが全て連続するようにもしたい*/

        /*横方向の座標が0の部分のアドレスを格納するメモリ領域を確保*/
        data = new RGB*[height];

        /*色データ用のメモリ領域全体を(アドレスが連続になるように)ここで確保して先頭アドレスをdata[0]に書き込む*/
        data[0] = new RGB[width*height];
        
        /*data[1], data[2], data[3],...,data[height_ - 1]にアドレスを書き込む*/
        for(int j=1; j<height; j++)
        {
            data[j] /*==&(data[j][0])*/ = data[0] + width * j;
        }
    }
    Image(): Image(input<int>("縦幅を入力してください。"), input<int>("横幅を入力してください。")){}


    /*色情報を書き込む*/ 
    constexpr void setpixel(int const& x, int const& y, RGB const& color) const
    {
        /*範囲内なら書き込む*/
        if((0<=x)&&(x<width)&&(0<=y)&&(y<height))
            data[y][x] = color;
    }

    /*色情報を読み込む*/
    constexpr RGB getpixel(int const& x, int const& y) const
    {
        /*範囲内ならその色が戻り値で、範囲外なら黒*/
        if((0<=x)&&(x<width)&&(0<=y)&&(y<height))
            return data[y][x];
        else
            return {0,0,0};
    }


    ~Image()
    {   
        // data→*dataの順にアクセスされるので、
        //dataを先に解放すると*dataを解放できなくなる

        /*色データ用領域を解放*/
        delete [] *data;

        /*アドレス用領域を解放*/
        delete [] data;
    }

    /*PNGファイルを作成(引数にファイル名)*/
    void output(std::string const& filename) const
    {
        writePNG(filename, width, height, *data);
    }

    /*PNGファイルを作成(標準入力からファイル名)*/
    void output() const
    {
        std::string filename = input<std::string>("ファイル名を入力してください。");
        writePNG(filename, width, height, *data);
    }
};



class Rainbow
{
    private:
    /*虹色を保存する配列*/
    short static constexpr color_num = 256;
    static RGB color[color_num];

    /*ベルカーブの分布でRGBを決める*/
    static constexpr unsigned char top_r = 206;
    static constexpr unsigned char top_g = 128;
    static constexpr unsigned char top_b = 50;

    /*ベルカーブの幅*/
    static constexpr unsigned int width = 1200;

    constexpr void cal_color() const
    {
        #define bell(s, a, x) (exp(-(double)(x-a)*(x-a)/(2*s)))
        for(int i=0; i<color_num; i++)
        {
            color[i] = {(char)(255*bell(width, top_r, i)), (char)(255*bell(width, top_g, i)), (char)(255*bell(width, top_b, i))};
        }
    }

    public:

    constexpr Rainbow()
    {
        cal_color();
    }

    RGB operator[](int const& i) const
    {
        if((0<=i)&&(i<color_num))
            return color[i];
        else if(i<0)
            return color[0];
        else
            return color[color_num - 1];
    }
};
RGB Rainbow::color[Rainbow::color_num];


////////////////////////////////////////////////////////////////
/*抽象クラス*/  //画素を書き込む処理を統一(子クラスでは座標から色を求める関数を実装するだけでいい)
class Pattern
{
    protected:
    Image& image;
    Pattern(Image& image_): image(image_){}

    public:
    /*Imageクラスに色情報を書き込むを子クラスで実装*/
    virtual RGB whatcolor(int const&, int const&) const = 0; 
    virtual void changesetting() = 0;

    constexpr void writecolor() const
    {
        for(int j = 0; j<image.height; j++)
        for(int i = 0; i<image.width; i++)
            image.setpixel(i, j, whatcolor(i, j));
    }
};

/*Patternから継承*/  //画像上の座標(i, j)から一様な色を求める処理を実装 
class Uniform : public Pattern
{
    public:
    RGB const color;

    Uniform(Image& image_, RGB color_) : Pattern(image_), color(color_){}
    Uniform(Image& image_): Uniform(image_, input<RGB>("RGBを空白区切りで入力してください。(0~255)")){}
    
    RGB whatcolor(int const& i, int const& j) const
    {
        return color;
    }

    void changesetting()
    {
        //何もしない
    }
};

/*Patternから継承*/  //画像上の座標(i, j)から縞模様を求める処理を実装 
class Sima : public Pattern
{
    public:
    int const simahaba;

    Sima(Image& image_, int simahaba_) : Pattern(image_), simahaba(simahaba_){}
    Sima(Image& image_): Sima(image_, input<int>("縞幅を入力してください。")){}

    RGB whatcolor(int const& i, int const& j) const
    {
        /*割り算の余りで縞模様を計算*/
        if((i/simahaba)%2 == 0)   
            return {(char)255,(char)255,(char)255};
        else
            return {0,0,0};
    }

    void changesetting()
    {
        //何もしない
    }
};


/*Patternから継承*/  //画像上の座標(i, j)から波の模様を求める処理を実装 
class Wave : public Pattern
{
    public:
    RGB const color;         //波の色
    double const wavelength; //波長
    Vec2 const direction;    //波の方向を向いた単位ベクトル

    Wave(Image& image_, RGB color_, int wavelength_, Vec2 direction_): Pattern(image_), color(color_), wavelength(wavelength_), direction(direction_.normalize()){}
    Wave(Image& image_): Wave(image_, input<RGB>("RGBを空白区切りで入力してください。(0~255)"), input<double>("波長を入力してください。"), input<Vec2>("波の方向を入力してください(空白区切りで)")){}

    RGB whatcolor(int const& i, int const& j) const
    {
        /*ある方向への距離を単位ベクトルとの内積で計算*/
        double distance = dot(direction, Vec2(i, j));

        /*波の位相*/
        double theta = (distance/wavelength) * 2 * M_PI;

        /*色の強さ(0~1)*/
        double r = (sin(theta) + 1)/2; //マイナスにならないように注意

        /*計算される色*/
        return color * r;
    }

    void changesetting()
    {
        //何もしない
    }
};

/*Patternから継承*/  //画像上の座標(i, j)からルールNの色を求める処理を実装 
class RuleN : public Pattern
{
    private:

    /*trueのところは黒, falseのところは白*/
    std::shared_ptr <bool[]> bits;

    /*char型から8つのビットを得る処理*/
    std::vector <bool> char_to_bit(char c)
    {
        std::vector<bool> b(8);
        for(int i=0; i<8; i++)
        {
            b[i] = c&1;
            c >>= 1;
        }
        return b;
    }

    /*事前に白か黒か計算して求める*/
    void calc_color()
    {
        for(int i=0; i<image.width ; i++)
        {
            if(i==image.width/2)
            {
                bits[i] = true;
            }
            else
            {
                bits[i] = false;
            }
        }

        for(int j=1; j<image.height; j++)
        for(int i=0; i<image.width ; i++)
        {
            int index = j*image.width + i;

            /*a, b, cはそれぞれ左上, 真上, 右上が黒ならtrue、そうでなければfalse*/ /*また、画素がなくてもfalse*/
            bool a, b, c;
            
            if(i == 0)
            {
                a = 0;
                b = bits[index - image.width - 0];
                c = bits[index - image.width + 1];
            }
            else if(i == image.width - 1)
            {
                a = bits[index - image.width - 1];
                b = bits[index - image.width - 0];
                c = 0;
            }
            else
            {
                a = bits[index - image.width - 1];
                b = bits[index - image.width - 0];
                c = bits[index - image.width + 1];
            }
            bits[index] = rule[4*a+2*b+c];
        }
    }

    public:
    /*ルールをビット列として求める*/
    std::vector <bool> const rule;


    RuleN(Image& image_, char rulenum_): Pattern(image_),
    bits{new bool[image_.width*image_.height]}, 
    rule(char_to_bit(rulenum_))
    {
        calc_color();
    }
    RuleN(Image& image_): RuleN(image_, input<short>("ルール番号を入力してください。(0-255)")){}


    /*画像上の座標(i, j)の色は何なのか求める関数*/
    RGB whatcolor(int const& i, int const& j) const
    {
        /*bits配列でtrueのところは黒, falseのところは白*/
        if(bits[image.width*j + i])
            return {0,0,0};
        else
            return {(char)255,(char)255,(char)255};
    }

    void changesetting()
    {
        //何もしない
    }
};

/*Patternから継承*/  //画像上の座標(i, j)から円形波の模様を求める処理を実装 
class CircWave : public Pattern
{
    public:
    RGB const color;
    Vec2 const center;
    int const wavelength;
    
    CircWave(Image& image_, RGB color_, Vec2 center_, int wavelength_): Pattern(image_), color(color_), center(center_), wavelength(wavelength_){}
    CircWave(Image& image_): CircWave(image_, input<RGB>("RGBを空白区切りで入力してください。(0~255)"), input<Vec2>("円の中心座標を空白区切りで入力してください"), input<int>("波長を入力してください")){}

    RGB whatcolor(int const& i, int const& j) const
    {
        /*円の中心からの距離*/
        double r = (Vec2(i, j) - center).length();

        /*波の位相を計算*/  /*円の中心での位相を0とした*/
        double theta = r/wavelength * 2 * M_PI;

        /*色の強さを計算*/
        double a = (cos(theta) + 1)/2;

        /*色を計算*/
        return color * a;
    }

    void changesetting()
    {
        //何もしない
    }
};

/*Patternから継承したが、まだ抽象クラス*/  //画像上の座標(i, j)から極座標の座標値を求める処理を実装 
class Polar : public Pattern
{
    protected:

    /*x軸と動径方向のなす角(0から2piの範囲)と距離を求める */
    constexpr void get_angle_distance(double& r, double& theta, int const& i, int const& j) const
    {
        /*v: 動径方向のベクトル*/
        Vec2 v = Vec2(i-center.x, center.y-j);

         /*distance: 渦の中心からの距離*/
        r = v.length();

        if(r == 0)
            theta = 0;
        else
        {
            v.normalize();
            if(v.y>=0)
                theta = acos(v.x);
            else
                theta = 2*M_PI - acos(v.x);
        }
    }

    public:
    Vec2 const center;
    Polar(Image& image_, Vec2 const& center_): Pattern(image_), center(center_){}
    Polar(Image& image_): Polar(image_, input<Vec2>("中心座標を空白区切りで入力してください")){}

    void changesetting()
    {
        //何もしない
    }
};

/*Polarから継承*/  //極座標の座標値から渦模様を求める処理を実装 
class GuruGuru_normal : public Polar
{
    public:
    int const simahaba;

    GuruGuru_normal(Image& image_, Vec2 const& center_, int const& simahaba_): Polar(image_, center_), simahaba(simahaba_){}
    GuruGuru_normal(Image& image_): GuruGuru_normal(image_, input<Vec2>("中心座標を空白区切りで入力してください"), input<int>("縞幅を入力してください")){}
    

    RGB whatcolor(int const& i, int const& j) const
    {
        /*x軸と動径方向のなす角(ラジアン)を求める (0から2piの範囲で)*/
        double r, theta;
        get_angle_distance(r, theta, i, j);

        /*渦を生成する処理*/
        if((int)(r/simahaba + theta/M_PI)%2 == 0)
            return {(char)255, (char)255, (char)255};
        else
            return {0, 0, 0};
    }
};

/*Polarから継承*/  //極座標の座標値から渦模様を求める処理を実装
class GuruGuru_EX : public Polar
{
    public:
    int const simahaba;

    GuruGuru_EX(Image& image_, Vec2 const& center_, int simahaba_): Polar(image_, center_), simahaba(simahaba_){}
    GuruGuru_EX(Image& image_): GuruGuru_EX(image_, input<Vec2>("中心座標を空白区切りで入力してください"), input<int>("渦の幅を入力してください")){}    

    RGB whatcolor(int const& i, int const& j) const
    {
        double r, theta;

        get_angle_distance(r, theta, i, j);

        bool f1, f2;

        if((int)(r/simahaba+10-(theta/M_PI*5 ))%2 == 0)
            f1=1;
        else
            f1=0;

        if(int(2+theta/(2*M_PI/20)-2*sin(r/500*(2*M_PI)*2))%2 == 0)
            f2=1;
        else
            f2=0;

        if(f1&&f2)
            return{(char)255, (char)255, 0};
        else if(f1&&!f2)
            return {(char)0, (char)0, (char)255};
        else if(!f1&&f2)
            return {(char)128, (char)128, (char)128};
        else
            return {(char)255, (char)0, (char)0};
    }
};

class Polar_rainbow : public Polar
{
    public:

    Polar_rainbow(Image& image_, Vec2 const& center_): Polar(image_, center_){}
    Polar_rainbow(Image& image_): Polar_rainbow(image_, input<Vec2>("中心座標を空白区切りで入力してください")){}

    RGB whatcolor(int const& i, int const& j) const
    {
        double r, theta;
        get_angle_distance(r, theta, i, j);

        static Rainbow rainbow;
        return rainbow[(int)theta/2/M_PI*255];
    }
};

/*Patternから継承したが、まだ抽象クラス*/  //画像上の座標(i, j)から模様を求める処理を実装したが、初期値を設定する処理仮想関数とした
class MJ : public Pattern
{
    protected:
    double length;
    Vec2 centorpos;
    double r;
    int samples;
    MJ(Image& image_): Pattern(image_)
    {
        samples = 100;
        length = 2;
        centorpos = Vec2(0, 0);
        r = 1;
    }

    /*マンデルブロ集合とジュリア集合では初期値が異なるのでinitという仮想関数を定義*/
    virtual void init(Complex&, Complex&, int const&, int const&) const = 0;


    public:

    /*画素の座標を複素平面上の座標に変換する処理*/
    constexpr void getuv(double& u, double& v, int const& i, int const& j)const
    {
        u = centorpos.x  +  length*(2*i - image.width) / image.width;
        v = centorpos.y  -  length*(2*j - image.height) / image.width;
    }

    /*マンデルブロ集合やジュリア集合を拡大・縮小・平行移動させるための処理*/
    void changesetting()
    {
        std::cout << "中心位置: "  << std::setprecision(15) << centorpos << std::endl;
        std::cout << "画像の横幅: " << length << "の2倍" << std::endl;

        double k = input<double>("何倍したいか入力してください。");
        length /= k;

        double rx = input<double>("横幅の何倍移動したいか入力してください。");
        centorpos.x += 2 * length * rx;
        double ry = input<double>("縦幅の何倍移動したいか入力してください。");
        centorpos.y += 2.0 * image.height/image.width * length * ry;   

        r = input<double>("カラー表示の比率を何倍にしたいか入力してください。");

        samples = input<int>("サンプル数を入力してください。");
    }
 
    RGB whatcolor(int const& i, int const& j) const
    {
        /*マンデルブロ集合とジュリア集合では初期値が異なるのでinitという仮想関数を定義*/
        Complex z, c;
        init(z, c, i, j);//初期値を定める関数

        /*diverge: 発散すればtrue,    num: 発散が判定されるまでの回数*/
        int num;
        bool diverge = false;
        for (num = 0; num < samples; ++num) 
        {
            z = z*z + c;
            /*絶対値が2を超えると発散することが知られている*/
            if (abs(z) > 2) 
            {
                diverge = true;
                break;
            }
        }
 
        /*色を求める*/
        if(diverge)
        {
            static Rainbow rainbow;
            return rainbow[(int)(255*(r*num/(double)samples))];
        }
        else
            return {0,0,0};
    }
};

/*MJから継承*/  //初期値を設定する処理(マンデルブロ集合用)を実装した
class Mandel : public MJ
{
    public:
    Mandel(Image& image_): MJ(image_){}

    /*initという仮想関数の実体(マンデルブロ集合用)をここに書く*/
    void init(Complex& z_, Complex& c_, int const& i, int const& j) const
    {
        /*マンデルブロ集合ではzの初期値は0で、cは座標から求める*/
        z_ = Complex(0, 0);
        getuv(c_.re, c_.im, i, j);
    }
};

/*MJから継承*/  //初期値を設定する処理(ジュリア集合用)を実装した
class Julia : public MJ
{
    public:
    /*cはメンバ*/
    Complex const c;

    Julia(Image& image_, Complex const& c_): MJ(image_), c(c_){}
    Julia(Image& image_): Julia(image_, input<Complex>("複素数を空白区切りで入力してください。")){}
    

    /*initという仮想関数の実体(ジュリア集合用)をここに書く*/
    void init(Complex& z_, Complex& c_, int const& i, int const& j) const
    {
        /*ジュリア集合ではzの初期値は位置から求め、cは任意*/
        getuv(z_.re, z_.im, i, j);
        c_ = c;//メンバの値をc_に渡す
    }
};


////////////////////////////////////////////////////////////////
/*ABC(クラス名)を空白区切りで書く*/
#define CLASSLIST ABC(Uniform) ABC(Sima) ABC(Wave) ABC(RuleN) ABC(CircWave) ABC(GuruGuru_normal) ABC(GuruGuru_EX) ABC(Polar_rainbow) ABC(Mandel) ABC(Julia)

#define ABC(a) a,
typedef enum class enm1
{
    CLASSLIST NUMBER
} PATTERNS;
#undef ABC


#define ABC(a) #a,
static const std::string patterns[(int)PATTERNS::NUMBER] = {CLASSLIST};
#undef ABC

/*
継承関係

Pattern......................画素を書き込む処理を統一(子クラスでは座標から色を求める関数を実装するだけでいい)
    ┣━Uniform 
    ┣━Sima        
    ┣━Wave
    ┣━RuleN
    ┣━CircWave   
    ┃
    ┣━━━Polar.............距離と角度を求める処理を統一(子クラスでは距離と角度から色を求める関数を実装するだけでいい)
    ┃     ┣━GuruGuru_normal
    ┃     ┣━GuruGuru_EX
    ┃     ┗━Polar_rainbow
    ┃
    ┃
    ┗━━━MJ................マンデルブロ集合, ジュリア集合を拡大・縮小する処理と発散判定の処理を統一(子クラスでは初期値の設定を実装するだけでいい)
        ┣━Mandel
        ┗━Julia

*/

/*動的ポリモーフィズムでPatternクラスの継承クラスをswitch文で条件分岐*/
std::shared_ptr <Pattern> setpattern(Image& image)
{
    std::cout << "-------------------------" << std::endl;

    for(long unsigned int i=0; i<sizeof(patterns)/sizeof(std::string); i++)
        std::cout << std::right  << std::setw(20) << patterns[i] + ": " << i << std::endl;

    std::cout << "-------------------------" << std::endl;

    int n = input<int>("番号を入力してください。");


    switch (n)
    {
        #define ABC(a) \
            case (int)PATTERNS::a: \
            return std::make_shared <a> (a(image)); \
            break;
        CLASSLIST
        #undef ABC
        
        default:
            std::cout << "異常な値が入力されました。" << std::endl;
            exit((int)ERRORCODE::SETPATTERN);
            break;
    }
}

constexpr void superimpose(Image const& from, Image& to, int const& x_begin, int const& y_begin, double const& magnification)
{
    int x_end = x_begin + magnification * to.width;
    int y_end = y_begin + magnification * to.height;

    for(int j = y_begin; j<y_end; j++)
    for(int i = x_begin; i<x_end; i++)
    {
        int x_from = (i - x_begin)/magnification;
        int y_from = (j - y_begin)/magnification;
        to.setpixel(i, j, from.getpixel(x_from, y_from));
    }
}


int main()
{
    while(1)
    {
        bool which;

        int width = input<int>("横幅を入力してください。");
        int height = input<int>("縦幅を入力してください。");
        Image image(width, height);
        std::shared_ptr <Pattern> pat = setpattern(image);

        while(1)
        {  
            pat->writecolor();
            image.output();
            which = input<bool>("同じ種類の画像を出力し直す場合は1を、やめる場合は0を入力してください。");
            if(which == 0)
                break;
            pat->changesetting();
        }

        which = input<bool>("最初からやり直す場合は1を、プログラムを終了する場合は0を入力してください。");
        if(which == 0)
            break;
    }
}