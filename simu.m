%%解析結果用プログラム

clear;
close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%以下パケットレベルの宣言

z=1744;%[bit]              %G.711 1パケット当たり218[Byte]=1744[bit]
%B=1000000;                     %1[Mbps]  1x10^6
B=1500000;
%b0=100000;%[bps]     https://www.school.ctc-g.co.jp/columns/inst/inst50.html

%μsで統一
param.K=5;
alphak=1/352000; % α1.0倍 [1/μs]
%alphak=1/281600; % α0.8倍 [1/μs]
betak=1/650000; % b1.0倍   [1/μs]
alphah=1/352000;
%alphah=1/281600;
betah=1/650000;
alphag=1/352000;
%alphag=1/281600;
betag=1/650000;
T=16000;         %  [μs]
N_list=1:1:20;
param.mu=B/(z*1000000)
%k:緊急呼,h:被災地呼,g:被災地外呼

lamdas0k=(betak)/(T*(alphak+betak));%param パケット到着率(3.1)
ca2jouk=(1-(1-alphak*T)^2)/((T^2)*(alphak+betak)^2); %ok 18.0950　平方変動係数(分散/平均)(3.2)
ca3jouk=ca2jouk^(3/2);
ca4jouk=ca2jouk^2;
sk0k=(2*alphak*T*(((alphak*T)^2)-3*alphak*T+3))/((alphak*T*(2-alphak*T))^(3/2)); %ひずみ度(3次中心積率/分散^3/2)(3.3)
Dk=(3*lamdas0k*(ca2jouk-1))/(2*sk0k*ca3jouk-3*ca4jouk-1);%param (3.6)
Fk=Dk*(3*ca4jouk-sk0k*ca3jouk-3*ca2jouk+2)/(3*(ca2jouk-1));%param (3.8)
Ek=Dk*(ca2jouk-1)/(Fk^2); %param (3.7)

lamdas0h=(betah)/(T*(alphah+betah)) ;
ca2jouh=(1-(1-alphah*T)^2)/((T^2)*(alphah+betah)^2); %ok 18.0950
ca3jouh=ca2jouh^(3/2);
ca4jouh=ca2jouh^2;
sk0h=(2*alphah*T*(((alphah*T)^2)-3*alphah*T+3))/((alphah*T*(2-alphah*T))^(3/2));
Dh=(3*lamdas0h*(ca2jouh-1))/(2*sk0h*ca3jouh-3*ca4jouh-1);
Fh=Dh*(3*ca4jouh-sk0h*ca3jouh-3*ca2jouh+2)/(3*(ca2jouh-1));
Eh=Dh*(ca2jouh-1)/(Fh^2);

lamdas0g=(betag)/(T*(alphag+betag));
ca2joug=(1-(1-alphag*T)^2)/((T^2)*(alphag+betag)^2); %ok 18.0950
ca3joug=ca2joug^(3/2);
ca4joug=ca2joug^2;
sk0g=(2*alphag*T*(((alphag*T)^2)-3*alphag*T+3))/((alphag*T*(2-alphag*T))^(3/2));
Dg=(3*lamdas0g*(ca2joug-1))/(2*sk0g*ca3joug-3*ca4joug-1);
Fg=Dg*(3*ca4joug-sk0g*ca3joug-3*ca2joug+2)/(3*(ca2joug-1));
Eg=Dg*(ca2joug-1)/(Fg^2);

%宣言
param.lamdas0k=lamdas0k
param.lamdas0h=lamdas0h
param.lamdas0g=lamdas0g
param.Dk=Dk;
param.Ek=Ek;
param.Fk=Fk;
param.Dh=Dh;
param.Eh=Eh;
param.Fh=Fh;
param.Dg=Dg;
param.Eg=Eg;
param.Fg=Fg;
param.Nk=0; 
param.Nh=0;
param.Ng=0;

tic; %ストップウォッチタイマーを開始

restmp=[]; 
res=0;
N_res=0;
ave=[];
num=0;
cou=0;
itmp=[];
jtmp=[];
ktmp=[];
high=[];
low=[];
ftmp=[];
flg=0;

res_ave=0;
res_temp=[];
prob_temp=[];
N=0;
kk=0;
mm=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%以下フローレベルの宣言

param.B    = 3;    %全帯域
param.b1   = 1;      %緊急呼占有帯域
param.b2   = 1;      %被災地呼占有帯域
param.b3   = 1;      %被災地外呼占有帯域

param.mu1  = 0.01; %緊急呼退去率
param.mu2  = 0.01; %被災地呼退去率
param.mu3  = 0.01; %被災地外呼退去率

c1=0.15; %緊急VoIPセッション上限呼損率
c2=0.5;  %被災地VoIPセッション上限呼損率

rho2d_list = [0.05:0.05:1];    %被災地VoIPセッションのトラヒック密度0.05から１まで0.05間隔でトラヒック密度をプロット
rho1d_list = [0.2];  %緊急VoIPセッションのトラヒック密度
rho3d_list = [0.2];  %被災地外VoIPセッションのトラヒック密度

Th2_list    = [0:param.B];  %閾値2を0から全帯域まで調べる thgin
Th1_list    = [0:param.B];  %閾値1を0から全帯域まで調べる thgout

Th1_opt_temp=[];
Th2_opt_temp=[];
opt_e=[];
opt_gin=[];
opt_gout=[];

H=1; %どういう役割？

thg=0;
thh=0;

tic; %時間測定


%原因は不明だが，固有値計算の途中でエラーになることがあるため，もし，エラーで止まった場合には，下記をコメントアウトして，
start_i =1;

for i=start_i:numel(rho1d_list)  %緊急呼のトラヒック密度
    for j=1:numel(rho2d_list)   %被災地呼のトラヒック密度
        for l=1:numel(rho3d_list)  %被災地外呼のトラヒック密度 
            %for x=H:numel(c2) %上限呼損率を求める場合に使用
                %c3=c2(H);
                th_gin=[];
                th_gout=[];
                e_call_late=[];
                gin_call_late=[];
                gout_call_late=[];
                packetloss=[];
                res_temp=[];
                prob_temp=[];
                a=1;
                b=1;
                min_gout=1;
                for k=1:numel(Th2_list)  %閾値2を1から全帯域まで 被災地外呼閾値thgin
                    a=k;
                    thh=k-1;
                    for m=1:numel(Th1_list)   %閾値1を1から(閾値2)-1まで　被災地呼閾値thgout
                        b=m;   
                        thg=m-1;
                        param.rho1dash = rho1d_list(i);
                        param.rho2dash = rho2d_list(j);
                        param.rho3dash = rho3d_list(l);
                        param.B_th2     = Th2_list(k);
                        param.B_th1     = Th1_list(m);
                        
                        [res_temp_e(k,m), res_temp_gin(k,m),res_temp_gout(k,m),totalpacketloss_temp(k,m)] = callbrock(param); %収容数後とのパケットロス率と各呼損率を取得
                                                
                        %各呼損率を各変数に格納
                        e_call_rate(k,m)=res_temp_e(k,m);
                        gin_call_rate(k,m)=res_temp_gin(k,m);
                        gout_call_rate(k,m)=res_temp_gout(k,m);                       
                        packetloss(k,m)=totalpacketloss_temp(k,m);  %パケットロスの期待値を配列として格納 

                        fprintf('th_gin:%d    th_gout:%d    res_ave=%f   r_e=%f    r_gin=%f    r_gout=%f\n',k,m,packetloss(k,m),e_call_rate(k,m),gin_call_rate(k,m),gout_call_rate(k,m));
                        

                        %被災地呼の呼損率が被災地外呼より小さく、被災地呼の呼損率がc２以下であり、緊急呼の呼損率がc１以下の時は各数値を変数に格納。そうでないときは０にする。
                        %%%目的関数のs.t.の条件
                        if(res_temp_gin(k,m)<res_temp_gout(k,m)&&res_temp_e(k,m)<=c1&&res_temp_gin(k,m)<=c2&&packetloss(k,m)<=0.0025&&k>m) 
                        %if(res_temp_gin(k,m)<res_temp_gout(k,m)&&res_temp_e(k,m)<=c1&&res_temp_gin(k,m)<=c2&&k>m)%通信品質を考慮しない場合
                            gout_call_rate(k,m)=res_temp_gout(k,m);
                            e_call_late(a,b)=res_temp_e(k,m);
                            gin_call_late(a,b)=res_temp_gin(k,m);
                            gout_call_late(a,b)=res_temp_gout(k,m);
 
                            th_gin(b)=m-1;
                            th_gout(a)=k-1; 

                        else
                                  
                            gout_call_rate(k,m) =0;
                            gout_call_late(a,b)=1;
                                                
                        end
                    end
                end
                
                
                
                if(res_temp_e(k,m)<c1&&res_temp_gin(k,m)<c2)           %緊急呼,被災地呼がある一定の呼損率以下なら
                    res_gout(k,m) =res_temp_gout(k,m);  %ある一定の呼損率以下の緊急呼の呼損率の同じ閾値での一般呼の呼損率
                else
                      res_gout(k,m)=1;           %それ以外は反映されないように
                end

                fprintf('rho1=%0.3f  rho2=%0.3f  rho3=%0.3f   \n', param.rho1dash, param.rho2dash, param.rho3dash);
                 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 最適閾値を抽出し、その時の呼損率を格納
                min_gout=min(gout_call_late,[],"all"); %%被災地外VoIPセッションの最小呼損率を抽出

                if(min_gout==1)  %最適閾値がない場合
                    Th2_opt_temp=[Th2_opt_temp nan];
                    Th1_opt_temp=[Th1_opt_temp nan];
                 
                    opt_e = [opt_e nan];  
                    opt_gin = [opt_gin nan];
                    opt_gout = [opt_gout nan];
                else    
                    [Th2_opt, Th1_opt] = find(gout_call_late == min_gout); %%その時の各閾値を抽出
                    Th2_opt_temp=[Th2_opt_temp Th2_opt-1];
                    Th1_opt_temp=[Th1_opt_temp Th1_opt-1];
                 
                    opt_e = [opt_e e_call_rate(Th2_opt,Th1_opt)];
                    opt_gin = [opt_gin gin_call_rate(Th2_opt,Th1_opt)];
                    opt_gout = [opt_gout gout_call_rate(Th2_opt,Th1_opt)];
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %end
                

                %被災地外VoIPセッション呼損率特性(x軸：被災地外VoIPセッション閾値、y軸：被災地VoIPセッション閾値、z軸：被災地外VoIPセッション呼損率) 
                 figure
                 x=Th1_list;  %%gout
                 y=Th2_list;  %%gin
                 z=gout_call_rate;
                 stem3(x,y,z,'r','o','LineStyle','none')
                 zlim([0 1.0])
                 xlabel('th_{gout}','FontSize',18)
                 ylabel('th_{gin}','FontSize',18)
                 zlabel('Generalout call blocking late,r_{gout}','FontSize',18)

                %パケットロス率特性(x軸：被災地外VoIPセッション閾値、y軸：被災地VoIPセッション閾値、z軸：パケットロス率の期待値) 
                 figure
                 for m=0:param.B  
                     mm=m+1; 
                     for k=0:param.B  
                          kk=k+1;
                          if(mm>=kk)
                              if(packetloss(mm,kk)>=0.0025) %パケットロス率が規定値以上であれば(通信品質を保証できていなければ)
                                s=scatter3(k,m,packetloss(mm,kk),48,[1 0 0],'LineWidth',2.0);  %%red 
                                s.SizeData=100;
                                hold on
                              else
                                s=scatter3(k,m,packetloss(mm,kk),48,[0 0 1],'LineWidth',2.0);  %%blue
                                s.SizeData=100;
                                hold on
                              end
                          end
                     end
                 end

                 hold off
                 xlabel('th_{gout}','FontSize',18)
                 ylabel('th_{gin}','FontSize',18)
                 zlabel('PacketLoss,L','FontSize',18)
      
        end
     
    end
    % iのループの最後で毎回ワークスペース中の全変数を保存
    save temp.mat

     

end


%%



toc;


 

 %最適閾値呼損率特性(横軸：被災地VoIPセッショントラヒック密度、縦軸：最適閾値のときの各呼損率)
 figure(1);
 plot(rho2d_list,opt_e, 'rp',rho2d_list,opt_gin, 'go',rho2d_list,opt_gout, 'b*','MarkerSize',15);
 xlabel('Traffic intencity \rho_2','FontSize',18)
 ylabel('call block probability','FontSize',18)
 legend({'e','gin','gout'},'FontSize',15)
 xlim([0, 1])
 ylim([0, 1])
 xticks(0:0.05:1)
 
 
 %最適閾値(横軸：被災地VoIPセッショントラヒック密度、縦軸：各最適閾値)
 figure(2);
 plot(rho2d_list,Th2_opt_temp, 'go',rho2d_list,Th1_opt_temp, 'b*','MarkerSize',15);
 xlabel('Traffic intencity \rho_2','FontSize',18)
 ylabel('optimal threshold th_{opt}','FontSize',18)
 legend({'thgin','thgout'},'FontSize',15)
 xlim([0, 1])
 ylim([0, param.B])
 xticks(0:0.05:1)
 yticks(0:1:param.B)
 