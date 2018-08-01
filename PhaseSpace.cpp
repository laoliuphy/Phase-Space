//
//  PhaseSpace.cpp
//  bb
//
//  Created by Yandong Liu on 3/14/16.
//  Copyright © 2016 Yandong Liu. All rights reserved.
//

#include "PhaseSpace.hpp"

/*two body phase space: two free parameters: theta and phi (x1, x2) theta=PI*x1 phi=2*PI*x2
 vector<double> b transfer the randam number
 vector<Particle> a save the particle four momentum information in the cm frame
 while the momentum information in the lab frame is saved in c*/
void Two_Body_P(vector<Particle>& a, vector<Particle>& c , vector<double>& b){
    Particle p1,p2, p, k1, k2;
    Particle bp1, bp2, bk1, bk2;
    p1=a[0];
    p2=a[1];
    k1=a[2];
    k2=a[3];
    p.SetPxPyPzE(0, 0, 0, 0);
    p+=p1;
    p+=p2;
    //   cout<<p.M()<<endl;
    double Q = sqrt(abs(lambda((p1+p2).M(), k1.M(), k2.M())))/(2*p.M());
    //   cout<<Q<<endl;
    double theta = PI_*b[0];
    double phi = 2*PI_*b[1];
    k1.SetPxPyPzE(Q*sin(theta)*cos(phi), Q*sin(theta)*sin(phi), Q*cos(theta), sqrt(Q*Q+k1.M2()));
    k2.SetPxPyPzE(Q*sin(PI_-theta)*cos(PI_+phi), Q*sin(PI_-theta)*sin(PI_+phi), Q*cos(PI_-theta), sqrt(Q*Q+k2.M2()));
    bk1=k1;
    bk2=k2;
    bp1=p1;
    bp2=p2;
    TVector3 boost;
    boost=p.BoostVector();
    bk1.Boost(boost);
    bk2.Boost(boost);
    
    c[0]=p1;
    c[1]=p2;
    c[2]=bk1;
    c[3]=bk2;//library frame
    
    a[2]=k1;
    a[3]=k2;
    //boot bp1 bp2 to cm frame
    bp1.Boost(-boost);
    bp2.Boost(-boost);
    a[0]=bp1;
    a[1]=bp2;
    
}

void Two_Body_P(vector<Particle>& a, vector<Particle>& c , vector<double>& b , double m , double width){
    Particle p1,p2, p, k1, k2;
    Particle bp1, bp2, bk1, bk2;
    p1=a[0];
    p2=a[1];
    k1=a[2];
    k2=a[3];
    p.SetPxPyPzE(0, 0, 0, 0);
    p+=p1;
    p+=p2;
    //   cout<<p.M()<<endl;
    double Q = sqrt(abs(lambda((p1+p2).M(), k1.M(), k2.M())))/(2*p.M());
    //   cout<<Q<<endl;
    double theta = PI_*b[0];
    double phi = 2*PI_*b[1];
    k1.SetPxPyPzE(Q*sin(theta)*cos(phi), Q*sin(theta)*sin(phi), Q*cos(theta), sqrt(Q*Q+k1.M2()));
    k2.SetPxPyPzE(Q*sin(PI_-theta)*cos(PI_+phi), Q*sin(PI_-theta)*sin(PI_+phi), Q*cos(PI_-theta), sqrt(Q*Q+k2.M2()));
    bk1=k1;
    bk2=k2;
    bp1=p1;
    bp2=p2;
    TVector3 boost;
    boost=p.BoostVector();
    bk1.Boost(boost);
    bk2.Boost(boost);
    
    c[0]=p1;
    c[1]=p2;
    c[2]=bk1;
    c[3]=bk2;//library frame
    
    a[2]=k1;
    a[3]=k2;
    //boot bp1 bp2 to cm frame
    bp1.Boost(-boost);
    bp2.Boost(-boost);
    a[0]=bp1;
    a[1]=bp2;
    
}

void Two_Body_P(int a, vector<Particle>& cm, vector<Particle>& lab, vector<double>& b){


    if (a==1) {
        if (cm[0].Mag()<10e-10) {
            Particle p, k1, k2;
            k1=cm[1];
            k2=cm[2];
            p=cm[0];
            double Q = sqrt(abs(lambda(p.M(), k1.M(), k2.M())))/(2*p.M());
            //   cout<<Q<<endl;
            double theta = PI_*b[0];
            double phi = 2*PI_*b[1];
            k1.SetPxPyPzE(Q*sin(theta)*cos(phi), Q*sin(theta)*sin(phi), Q*cos(theta), sqrt(Q*Q+k1.M2()));
            k2.SetPxPyPzE(Q*sin(PI_-theta)*cos(PI_+phi), Q*sin(PI_-theta)*sin(PI_+phi), Q*cos(PI_-theta), sqrt(Q*Q+k2.M2()));
            cm[0]=p;
            cm[1]=k1;
            cm[2]=k2;
            lab[0]=p;
            lab[1]=k1;
            lab[2]=k2;
        } else {
            Particle p, k1, k2;
            Particle labp, labk1, labk2;
            
            p=cm[0];
            labp=cm[0];
            k1=cm[1];
            k2=cm[2];
            double Q = sqrt(abs(lambda(labp.M(), k1.M(), k2.M())))/(2*labp.M());
            //   cout<<Q<<endl;
            double theta = PI_*b[0];
            double phi = 2*PI_*b[1];
            k1.SetPxPyPzE(Q*sin(theta)*cos(phi), Q*sin(theta)*sin(phi), Q*cos(theta), sqrt(Q*Q+k1.M2()));
            k2.SetPxPyPzE(Q*sin(PI_-theta)*cos(PI_+phi), Q*sin(PI_-theta)*sin(PI_+phi), Q*cos(PI_-theta), sqrt(Q*Q+k2.M2()));
            labk1=k1;
            labk2=k2;
            TVector3 boost;
            boost=p.BoostVector();
            p.Boost(-boost);
            labk1.Boost(boost);
            labk2.Boost(boost);
            cm[0]=p;
            cm[1]=k1;
            cm[2]=k2;
            lab[0]=p;
            lab[1]=k1;
            lab[2]=k2;
        }
    } else if (a==2){
        Two_Body_P(cm, lab, b);
    }
    
}


/*p1+p2->k1+q1->k1+k2+k3
 5 free parameters k2 k3 system theta2 phi2
 k1 q1 system theta1 phi1 and the effective mass of q1 */
void Three_Body_P(vector<Particle>& a, vector<Particle>& c, vector<double>& b){
    Particle p1,p2, k1, k2, k3, q1;
    Particle qk1, qk2, qk3, qp1, qp2;
    p1=a[0];
    p2=a[1];
    k1=a[2];
    k2=a[3];
    k3=a[4];
    qp1=p1;
    qp2=p2;
    double tmp_m = k2.M()+k3.M() + ((p1+p2).M()-k1.M()-k2.M()-k3.M())*b[0];
    double theta1 = b[1]*PI_;
    double phi1 = b[2]*2*PI_;
    double theta2 = b[3]*PI_;
    double phi2 = b[4]*2*PI_;
    double Q1 = sqrt(abs(lambda((p1+p2).M(), k1.M(), tmp_m)))/(2*(p1+p2).M()); // momentum of k1 and q1 systme
    double Q2 = sqrt(abs(lambda(tmp_m, k2.M(), k3.M())))/(2*tmp_m); // momentum of k2 and k3 system
    k1.SetPxPyPzE(Q1*sin(theta1)*cos(phi1), Q1*sin(theta1)*sin(phi1), Q1*cos(theta1), sqrt(Q1*Q1+k1.M2()));
    k2.SetPxPyPzE(Q2*sin(theta2)*cos(phi2), Q2*sin(theta2)*sin(phi2), Q2*cos(theta2), sqrt(Q2*Q2+k2.M2()));
    k3.SetPxPyPzE(Q2*sin(PI_-theta2)*cos(PI_+phi2), Q2*sin(PI_-theta2)*sin(PI_+phi2), Q2*cos(PI_-theta2), sqrt(Q2*Q2+k3.M2()));
    q1.SetPxPyPzE(Q1*sin(PI_-theta1)*cos(PI_+phi1), Q1*sin(PI_-theta1)*sin(PI_+phi1), Q1*cos(PI_-theta1), sqrt(Q1*Q1+tmp_m*tmp_m));
    
    TVector3 boost;
    boost=q1.BoostVector();
    
    k2.Boost(boost);
    k3.Boost(boost); // boost the k2 k3 to the q1 kinematic frame (p1 p2 cm frame)
    
    qk1=k1;
    qk2=k2;
    qk3=k3;
    
    TVector3 boostp;
    boostp=(p1+p2).BoostVector();
    
    qk1.Boost(boostp);
    qk2.Boost(boostp);
    qk3.Boost(boostp); //  boost the k1 k2 k3 to the laboratory frame
    
    qp1.Boost(-boostp);
    qp2.Boost(-boostp); // boost to cm frame
    
    a[0]=qp1;
    a[1]=qp2;
    a[2]=k1;
    a[3]=k2;
    a[4]=k3; // cm frame
    
    c[0]=p1;
    c[1]=p2;
    c[2]=qk1;
    c[3]=qk2;
    c[4]=qk3; //laboratory frame
//    cout<<"***************"<<endl;
//    cout<<tmp_m<<" "<<(a[3]+a[4]).M()<<" "<<(c[3]+c[4]).M()<<endl;
//    cout<<"***************"<<endl;
}

void Three_Body_P(vector<Particle>& a, vector<Particle>& c , vector<double>& b, double mass, double width){

    if ((a[0]+a[1]).M()<mass) {
        Three_Body_P(a, c, b);
    } else {
        Particle p1,p2, k1, k2, k3, q1;
        Particle qk1, qk2, qk3, qp1, qp2;
        p1=a[0];
        p2=a[1];
        k1=a[2];
        k2=a[3];
        k3=a[4];
        qp1=p1;
        qp2=p2;
        double lower=k2.M()+k3.M();
        double higher=(p1+p2).M()-k1.M();
        double tlower=atan((lower*lower-mass*mass)/(mass*width));
        double thigher=atan((higher*higher-mass*mass)/(mass*width));
        double tmass=tlower+(thigher-tlower)*b[0];
        double tmp_m=sqrt(mass*width*tan(tmass)+mass*mass);
        double theta1 = b[1]*PI_;
        double phi1 = b[2]*2*PI_;
        double theta2 = b[3]*PI_;
        double phi2 = b[4]*2*PI_;
        double Q1 = sqrt(abs(lambda((p1+p2).M(), k1.M(), tmp_m)))/(2*(p1+p2).M()); // momentum of k1 and q1 systme
        double Q2 = sqrt(abs(lambda(tmp_m, k2.M(), k3.M())))/(2*tmp_m); // momentum of k2 and k3 system
        k1.SetPxPyPzE(Q1*sin(theta1)*cos(phi1), Q1*sin(theta1)*sin(phi1), Q1*cos(theta1), sqrt(Q1*Q1+k1.M2()));
        k2.SetPxPyPzE(Q2*sin(theta2)*cos(phi2), Q2*sin(theta2)*sin(phi2), Q2*cos(theta2), sqrt(Q2*Q2+k2.M2()));
        k3.SetPxPyPzE(Q2*sin(PI_-theta2)*cos(PI_+phi2), Q2*sin(PI_-theta2)*sin(PI_+phi2), Q2*cos(PI_-theta2), sqrt(Q2*Q2+k3.M2()));
        q1.SetPxPyPzE(Q1*sin(PI_-theta1)*cos(PI_+phi1), Q1*sin(PI_-theta1)*sin(PI_+phi1), Q1*cos(PI_-theta1), sqrt(Q1*Q1+tmp_m*tmp_m));
        
        TVector3 boost;
        boost=q1.BoostVector();
        
        k2.Boost(boost);
        k3.Boost(boost); // boost the k2 k3 to the q1 kinematic frame (p1 p2 cm frame)
        
        qk1=k1;
        qk2=k2;
        qk3=k3;
        
        TVector3 boostp;
        boostp=(p1+p2).BoostVector();
        
        qk1.Boost(boostp);
        qk2.Boost(boostp);
        qk3.Boost(boostp); //  boost the k1 k2 k3 to the laboratory frame
        
        qp1.Boost(-boostp);
        qp2.Boost(-boostp); // boost to cm frame
        
        a[0]=qp1;
        a[1]=qp2;
        a[2]=k1;
        a[3]=k2;
        a[4]=k3; // cm frame
        
        c[0]=p1;
        c[1]=p2;
        c[2]=qk1;
        c[3]=qk2;
        c[4]=qk3; //laboratory frame
    }
}

void Three_Body_P(int n, vector<Particle>& cm, vector<Particle>& lab, vector<double>& b){
    if (n==1) {
        if (cm[0].Mag()<10e-10) {
            Particle p, k1, k2, k3, q1;
            p=cm[0];
            k1=cm[1];
            k2=cm[2];
            k3=cm[3];
            double tmp_m = k2.M()+k3.M() + (p.M()-k1.M()-k2.M()-k3.M())*b[0];
            double theta1 = b[1]*PI_;
            double phi1 = b[2]*2*PI_;
            double theta2 = b[3]*PI_;
            double phi2 = b[4]*2*PI_;
            double Q1 = sqrt(abs(lambda(p.M(), k1.M(), tmp_m)))/(2*p.M()); // momentum of k1 and q1 systme
            double Q2 = sqrt(abs(lambda(tmp_m, k2.M(), k3.M())))/(2*tmp_m); // momentum of k2 and k3 system
            k1.SetPxPyPzE(Q1*sin(theta1)*cos(phi1), Q1*sin(theta1)*sin(phi1), Q1*cos(theta1), sqrt(Q1*Q1+k1.M2()));
            k2.SetPxPyPzE(Q2*sin(theta2)*cos(phi2), Q2*sin(theta2)*sin(phi2), Q2*cos(theta2), sqrt(Q2*Q2+k2.M2()));
            k3.SetPxPyPzE(Q2*sin(PI_-theta2)*cos(PI_+phi2), Q2*sin(PI_-theta2)*sin(PI_+phi2), Q2*cos(PI_-theta2), sqrt(Q2*Q2+k3.M2()));
            q1.SetPxPyPzE(Q1*sin(PI_-theta1)*cos(PI_+phi1), Q1*sin(PI_-theta1)*sin(PI_+phi1), Q1*cos(PI_-theta1), sqrt(Q1*Q1+tmp_m*tmp_m));
            TVector3 boost;
            boost=q1.BoostVector();
            k2.Boost(boost);
            k3.Boost(boost); // boost the k2 k3 to the q1 kinematic frame ( p frame)
            cm[0]=p;
            cm[1]=k1;
            cm[2]=k2;
            cm[3]=k3;
            lab=cm;
        } else {
            Particle p, k1, k2, k3, q1;
            Particle qk1, qk2, qk3, qp;
            p=cm[0];
            k1=cm[1];
            k2=cm[2];
            k3=cm[3];
            qp=cm[0];
            double tmp_m = k2.M()+k3.M() + (p.M()-k1.M()-k2.M()-k3.M())*b[0];
            double theta1 = b[1]*PI_;
            double phi1 = b[2]*2*PI_;
            double theta2 = b[3]*PI_;
            double phi2 = b[4]*2*PI_;
            double Q1 = sqrt(abs(lambda(p.M(), k1.M(), tmp_m)))/(2*p.M()); // momentum of k1 and q1 systme
            double Q2 = sqrt(abs(lambda(tmp_m, k2.M(), k3.M())))/(2*tmp_m); // momentum of k2 and k3 system
            k1.SetPxPyPzE(Q1*sin(theta1)*cos(phi1), Q1*sin(theta1)*sin(phi1), Q1*cos(theta1), sqrt(Q1*Q1+k1.M2()));
            k2.SetPxPyPzE(Q2*sin(theta2)*cos(phi2), Q2*sin(theta2)*sin(phi2), Q2*cos(theta2), sqrt(Q2*Q2+k2.M2()));
            k3.SetPxPyPzE(Q2*sin(PI_-theta2)*cos(PI_+phi2), Q2*sin(PI_-theta2)*sin(PI_+phi2), Q2*cos(PI_-theta2), sqrt(Q2*Q2+k3.M2()));
            q1.SetPxPyPzE(Q1*sin(PI_-theta1)*cos(PI_+phi1), Q1*sin(PI_-theta1)*sin(PI_+phi1), Q1*cos(PI_-theta1), sqrt(Q1*Q1+tmp_m*tmp_m));
            TVector3 boost;
            boost=q1.BoostVector();
            
            k2.Boost(boost);
            k3.Boost(boost); // boost the k2 k3 to the q1 kinematic frame (p1 p2 cm frame)
           
            qk1=k1;
            qk2=k2;
            qk3=k3;
            
            TVector3 boostp;
            boostp=p.BoostVector();
           
            qk1.Boost(boostp);
            qk2.Boost(boostp);
            qk3.Boost(boostp); //  boost the k1 k2 k3 to the laboratory frame
            
            lab[0]=p;
            lab[1]=qk1;
            lab[2]=qk2;
            lab[3]=qk3;
            
            qp.Boost(-boostp);
            cm[0]=qp;
            cm[1]=k1;
            cm[2]=k2;
            cm[3]=k3;
        }
    } else if (n==2){
        Three_Body_P(cm,lab, b);
    }
    
}// n is the number of the initial particle


void Three_Body_P(int n, vector<Particle>& cm, vector<Particle>& lab, vector<double>& b, double mass, double width){
    if (n==1) {
        if (cm[0].Mag()<10e-10) {
            Particle p, k1, k2, k3, q1;
            p=cm[0];
            k1=cm[1];
            k2=cm[2];
            k3=cm[3];
            double lower=k2.M()+k3.M();
            double higher=p.M()-k1.M();
            double tlower=atan((lower*lower-mass*mass)/(mass*width));
            double thigher=atan((higher*higher-mass*mass)/(mass*width));
            double tmass=tlower+(thigher-tlower)*b[0];
            double tmp_m=sqrt(mass*width*tan(tmass)+mass*mass);
            double theta1 = b[1]*PI_;
            double phi1 = b[2]*2*PI_;
            double theta2 = b[3]*PI_;
            double phi2 = b[4]*2*PI_;
            double Q1 = sqrt(abs(lambda(p.M(), k1.M(), tmp_m)))/(2*p.M()); // momentum of k1 and q1 systme
            double Q2 = sqrt(abs(lambda(tmp_m, k2.M(), k3.M())))/(2*tmp_m); // momentum of k2 and k3 system
            k1.SetPxPyPzE(Q1*sin(theta1)*cos(phi1), Q1*sin(theta1)*sin(phi1), Q1*cos(theta1), sqrt(Q1*Q1+k1.M2()));
            k2.SetPxPyPzE(Q2*sin(theta2)*cos(phi2), Q2*sin(theta2)*sin(phi2), Q2*cos(theta2), sqrt(Q2*Q2+k2.M2()));
            k3.SetPxPyPzE(Q2*sin(PI_-theta2)*cos(PI_+phi2), Q2*sin(PI_-theta2)*sin(PI_+phi2), Q2*cos(PI_-theta2), sqrt(Q2*Q2+k3.M2()));
            q1.SetPxPyPzE(Q1*sin(PI_-theta1)*cos(PI_+phi1), Q1*sin(PI_-theta1)*sin(PI_+phi1), Q1*cos(PI_-theta1), sqrt(Q1*Q1+tmp_m*tmp_m));
            TVector3 boost;
            boost=q1.BoostVector();
            k2.Boost(boost);
            k3.Boost(boost); // boost the k2 k3 to the q1 kinematic frame ( p frame)
            cm[0]=p;
            cm[1]=k1;
            cm[2]=k2;
            cm[3]=k3;
            lab=cm;
        } else {
//            cout<<"come here!"<<endl;
            Particle p, k1, k2, k3, q1;
            Particle qk1, qk2, qk3, qp;
            p=cm[0];
            k1=cm[1];
            k2=cm[2];
            k3=cm[3];
            qp=cm[0];
            double lower=k2.M()+k3.M();
            double higher=p.M()-k1.M();
            double tlower=atan((lower*lower-mass*mass)/(mass*width));
            double thigher=atan((higher*higher-mass*mass)/(mass*width));
            double tmass=tlower+(thigher-tlower)*b[0];
            double tmp_m=sqrt(mass*width*tan(tmass)+mass*mass);
//            cout<<tlower<<" "<<thigher<<" "<<tmp_m<<endl;
            double theta1 = b[1]*PI_;
            double phi1 = b[2]*2*PI_;
            double theta2 = b[3]*PI_;
            double phi2 = b[4]*2*PI_;
            double Q1 = sqrt(abs(lambda(p.M(), k1.M(), tmp_m)))/(2*p.M()); // momentum of k1 and q1 systme
            double Q2 = sqrt(abs(lambda(tmp_m, k2.M(), k3.M())))/(2*tmp_m); // momentum of k2 and k3 system
            k1.SetPxPyPzE(Q1*sin(theta1)*cos(phi1), Q1*sin(theta1)*sin(phi1), Q1*cos(theta1), sqrt(Q1*Q1+k1.M2()));
            k2.SetPxPyPzE(Q2*sin(theta2)*cos(phi2), Q2*sin(theta2)*sin(phi2), Q2*cos(theta2), sqrt(Q2*Q2+k2.M2()));
            k3.SetPxPyPzE(Q2*sin(PI_-theta2)*cos(PI_+phi2), Q2*sin(PI_-theta2)*sin(PI_+phi2), Q2*cos(PI_-theta2), sqrt(Q2*Q2+k3.M2()));
            q1.SetPxPyPzE(Q1*sin(PI_-theta1)*cos(PI_+phi1), Q1*sin(PI_-theta1)*sin(PI_+phi1), Q1*cos(PI_-theta1), sqrt(Q1*Q1+tmp_m*tmp_m));
            TVector3 boost;
            boost=q1.BoostVector();
            
            k2.Boost(boost);
            k3.Boost(boost); // boost the k2 k3 to the q1 kinematic frame (p1 p2 cm frame)
            
            qk1=k1;
            qk2=k2;
            qk3=k3;
            
            TVector3 boostp;
            boostp=p.BoostVector();
            
            qk1.Boost(boostp);
            qk2.Boost(boostp);
            qk3.Boost(boostp); //  boost the k1 k2 k3 to the laboratory frame
            
            lab[0]=p;
            lab[1]=qk1;
            lab[2]=qk2;
            lab[3]=qk3;
            
            qp.Boost(-boostp);
            cm[0]=qp;
            cm[1]=k1;
            cm[2]=k2;
            cm[3]=k3;
//            cout<<lab[0].M()<<endl;
        }
    } else if (n==2){
        Three_Body_P(cm, lab, b, mass, width);
    }
}


double Two_Body(vector<Particle>& a , vector<Particle>& c , vector<double>& b){
    
    Two_Body_P(a, c, b);
    
    //    cout<<a[0].Px()<<" "<<a[0].Py()<<" "<<a[0].Pz()<<" "<<a[0].E()<<endl;
    //    cout<<a[1].Px()<<" "<<a[1].Py()<<" "<<a[1].Pz()<<" "<<a[1].E()<<endl;
    //   cout<<(a[0]+a[1]).Px()<<" "<<(a[0]+a[1]).Py()<<" "<<(a[0]+a[1]).Pz()<<" "<<(a[0]+a[1]).E()<<endl;
    //    cout<<"~~~~~~~~~~~~~~"<<endl;
    double tmp = (0.125*sin(b[0]*PI_)*a[2].P())/(a[0]+a[1]).M();
    //cout<<a[2].P()<<endl;
    //cout<<(a[0]+a[1]).M()<<endl;
   // cout<<tmp<<endl;
    return tmp;
    
}

double Two_Body(int n, vector<Particle>& a, vector<Particle>& c, vector<double>& b){

    Two_Body_P(n, a, c, b);
    double tmp=0;
    if (n==1) {
        tmp=(0.125*sin(b[0]*PI_)*a[2].P())/a[0].M();
    } else if (n==2){
        tmp=(0.125*sin(b[0]*PI_)*a[2].P())/(a[0]+a[1]).M();
    }
    
    return tmp;
}

double Three_Body(vector<Particle> & cm , vector<Particle> & lab, vector<double> & b){
    Three_Body_P(cm, lab, b);
    Particle p,q;
    q=cm[3];
    q+=cm[4];
    p=cm[0];
    p+=cm[1];
    vector<Particle> phase1(3), phaselab1(3);
    vector<Particle> phase2(3), phaselab2(3);
    phase1[0]=p;
    phase1[1]=cm[2];
    phase1[2]=q;
    phase2[0]=q;
    phase2[1]=cm[3];
    phase2[2]=cm[4];
    
    
    
    double tmp=0, tmp1=0, tmp2=0, tmp3=0;
    
    vector<double> r1(2);
    vector<double> r2(2);
    r1[0]=b[1];
    r1[1]=b[2];
    r2[0]=b[3];
    r2[1]=b[4];
    
    tmp1=Two_Body(1, phase1, phaselab1 , r1);
    
//    cout<<"come here!"<<endl;
    tmp2=Two_Body(1, phase2, phaselab2, r2);
 
    tmp3=0.5*(2*q.M()*(p.M()-cm[2].M()-cm[3].M()-cm[4].M()))/PI_;
//    cout<<tmp1<<" "<<tmp2<<" "<<tmp3<<endl;
//    cout<<p.M()-cm[2].M()-cm[3].M()-cm[4].M()<<endl;
    tmp=tmp1*tmp2*tmp3;
    return tmp;
    
}

double Three_Body(int n, vector<Particle> & cm, vector<Particle> & lab, vector<double> & b){
    double tmp=0;
    if (n==1) {
        Three_Body_P(n,cm, lab, b);
        Particle p,q;
        q=cm[2];
        q+=cm[3];
        p=cm[0];
        vector<Particle> phase1(3), phaselab1(3);
        vector<Particle> phase2(3), phaselab2(3);
        phase1[0]=p;
        phase1[1]=cm[1];
        phase1[2]=q;
        phase2[0]=q;
        phase2[1]=cm[2];
        phase2[2]=cm[3];
        
        
        
        double tmp1=0, tmp2=0, tmp3=0;
        
        vector<double> r1(2);
        vector<double> r2(2);
        r1[0]=b[1];
        r1[1]=b[2];
        r2[0]=b[3];
        r2[1]=b[4];
        
        tmp1=Two_Body(1, phase1, phaselab1 , r1);
        
        //    cout<<"come here!"<<endl;
        tmp2=Two_Body(1, phase2, phaselab2, r2);
        
        tmp3=0.5*(2*q.M()*(p.M()-cm[1].M()-cm[2].M()-cm[3].M()))/PI_;
        
        tmp=tmp1*tmp2*tmp3;
        
    } else if (n==2){
        tmp=Three_Body(cm, lab, b);
    }
    return tmp;
}

double Three_Body(int n, vector<Particle>& cm, vector<Particle>& lab, vector<double>& b ,double mass, double width){
    double tmp=0;
    if (n==1) {
 //       cout<<"come here!"<<endl;
        Three_Body_P(1, cm, lab, b, mass, width);
        Particle p,q;
        q=cm[2];
        q+=cm[3];
        p=cm[0];
        vector<Particle> phase1(3), phaselab1(3);
        vector<Particle> phase2(3), phaselab2(3);
        phase1[0]=p;
        phase1[1]=cm[1];
        phase1[2]=q;
        phase2[0]=q;
        phase2[1]=cm[2];
        phase2[2]=cm[3];
        
        
        
        double tmp1=0, tmp2=0, tmp3=0;
        
        vector<double> r1(2);
        vector<double> r2(2);
        r1[0]=b[1];
        r1[1]=b[2];
        r2[0]=b[3];
        r2[1]=b[4];
        
        tmp1=Two_Body(1, phase1, phaselab1 , r1);
        
        //    cout<<"come here!"<<endl;
        tmp2=Two_Body(1, phase2, phaselab2, r2);
        
        double lower=cm[2].M()+cm[3].M();
        double higher=cm[0].M()-cm[1].M();
        double tlower=atan((lower*lower-mass*mass)/(mass*width));
        double thigher=atan((higher*higher-mass*mass)/(mass*width));
        double tmass=tlower+(thigher-tlower)*b[0];
//        double tmp_m=sqrt(mass*width*tan(tmass)+mass*mass);
        
        tmp3=0.5*(mass*width*(thigher-tlower))/(PI_*cos(tmass)*cos(tmass));
        
//       cout<<tmp1<<" "<<tmp2<<" "<<tmp3<<endl;
        
        tmp=tmp1*tmp2*tmp3;
        
    } else if (n==2){
        Three_Body_P(2, cm, lab, b, mass, width);
        Particle p,q;
        q=cm[3];
        q+=cm[4];
        p=cm[0];
        p+=cm[1];
        vector<Particle> phase1(3), phaselab1(3);
        vector<Particle> phase2(3), phaselab2(3);
        phase1[0]=p;
        phase1[1]=cm[2];
        phase1[2]=q;
        phase2[0]=q;
        phase2[1]=cm[3];
        phase2[2]=cm[4];
        
        
        
        double tmp1=0, tmp2=0, tmp3=0;
        
        vector<double> r1(2);
        vector<double> r2(2);
        r1[0]=b[1];
        r1[1]=b[2];
        r2[0]=b[3];
        r2[1]=b[4];
        
        tmp1=Two_Body(1, phase1, phaselab1 , r1);
        
        //    cout<<"come here!"<<endl;
        tmp2=Two_Body(1, phase2, phaselab2, r2);
        
        double lower=cm[3].M()+cm[4].M();
        double higher=p.M()-cm[1].M();
        double tlower=atan((lower*lower-mass*mass)/(mass*width));
        double thigher=atan((higher*higher-mass*mass)/(mass*width));
        double tmass=tlower+(thigher-tlower)*b[0];
        //        double tmp_m=sqrt(mass*width*tan(tmass)+mass*mass);
        
        tmp3=0.5*(mass*width*(thigher-tlower))/(PI_*cos(tmass)*cos(tmass));
        
//        cout<<tmp1<<" "<<tmp2<<" "<<tmp3<<endl;
        
        
        tmp=tmp1*tmp2*tmp3;
//        cout<<tmp<<endl;
    }
    
//    cout<<tmp<<endl;
    return tmp;
}





