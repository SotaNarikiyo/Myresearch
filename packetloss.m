% 収容数等を引数としてパケットロス率を出力する関数
% 各呼の組み合わせにおけるパケットロス率の期待値を出力する関数

function totalr = packetloss(param)


%z=1744;%[bit]              %G.711 1パケット当たり218[Byte]=1744[bit]
%B=100000000;%[bps]
%b0=100000;%[bps]     https://www.school.ctc-g.co.jp/columns/inst/inst50.html
%b1=100;             %通話あたりの音声帯域(正確には87.2kbpsだが、実環境では多少余裕を見た値を想定することが一般的)
%b2=100;

K=param.K;   %%%%
mu=param.mu; %%%%

Nk=param.Nk;
Nh=param.Nh;
Ng=param.Ng;

lamdas0k=param.lamdas0k;
lamdas0h=param.lamdas0h;
lamdas0g=param.lamdas0g;
Dk=param.Dk;
Ek=param.Ek;
Fk=param.Fk;
Dh=param.Dh;
Eh=param.Eh;
Fh=param.Fh;
Dg=param.Dg;
Eg=param.Eg;
Fg=param.Fg;


r0k=Dk*(1+(1/sqrt(1+Nk*lamdas0k*Ek))); %r0 位相の変化率
r1k=Dk*(1-(1/sqrt(1+Nk*lamdas0k*Ek))); %r1 位相の変化率
r0h=Dh*(1+(1/sqrt(1+Nh*lamdas0h*Eh))); %
r1h=Dh*(1-(1/sqrt(1+Nh*lamdas0h*Eh)));
r0g=Dg*(1+(1/sqrt(1+Ng*lamdas0g*Eg)));
r1g=Dg*(1-(1/sqrt(1+Ng*lamdas0g*Eg)));

lamdak0=Nk*lamdas0k+Fk+Fk*sqrt(1+Nk*lamdas0k*Ek);
lamdak1=Nk*lamdas0k+Fk-Fk*sqrt(1+Nk*lamdas0k*Ek);
lamdah0=Nh*lamdas0h+Fh+Fh*sqrt(1+Nh*lamdas0h*Eh);
lamdah1=Nh*lamdas0h+Fh-Fh*sqrt(1+Nh*lamdas0h*Eh);
lamdag0=Ng*lamdas0g+Fg+Fg*sqrt(1+Ng*lamdas0g*Eg);
lamdag1=Ng*lamdas0g+Fg-Fg*sqrt(1+Ng*lamdas0g*Eg);


%%smatrix(スパース行列(疎行列))
% 状態番号を格納する行列smatrixと，状態ごとの（各帯域の）収容済み呼数を格納するstateをつくる
smatrix = ones(K+3,11)*NaN;  %%%%%%%%%%%%%%%%%,K+3,K+3,5  状態数は最大呼数+1，さらに前後に番兵を設置するためもう+2する 8次元の１の配列
statenum=0;


%キューにパケットを最大限収容する
for i=0:K
    for q=0:7                %それが5以下なら・・・
        statenum=statenum+1;    %まずは状態番号をインクリメント
        smatrix(i+2,q+2)=statenum;     % 本数0の状態があるので+1，番兵の分さらに+1
        state(statenum).packet=i;     %状態番号がstatenumのときの各位相のときのパケット数を保存
        state(statenum).q=q;          
    end
end



% Smatrixへのaccessor
statemat = @(i,q) smatrix(i+2,q+2);  %%%%%%%%%%%%%%%% %%これ以降，直接Smatrixは参照しないこと！

% 隣接行列Aをゼロで初期化
A=zeros(statenum,statenum);
A=sparse(A);

loss000=[];
loss001=[];
loss010=[];
loss011=[];
loss100=[];
loss101=[];
loss110=[];
loss111=[];


% 隣接行列Aの中身を作成
for s=1:statenum    %通し番号sをスタート地点とする

    packet=state(s).packet;             %000  現在の通し番号におけるキューで待機しているパケット数
    q=state(s).q;
    
    if(q==0)  %%斜めの遷移を制御
        %---q000の位相でのパケット到着        
        t = statemat(packet+1,0);  %%%%%%%%%%%%%%            %iを一つ増やしてみたときの通し番号を行き先tとする．
        if(~isnan(t))                       %上のような状態が存在するかチェック．存在しないならt=NaNのはず
            A(t,s)=lamdak0+lamdah0+lamdag0;                 %状態がある場合：スタート地点sから行き先tに遷移する確率はλ0
        else                                %網内マックス呼損
            loss000 = [loss000 s];              %状態が無い場合：いま到着した狭帯域は呼損するので，                                        
        end                                  %その番号を配列loss1に控えておく（後で呼損率の計算に使う）．
    end                                 
    

    if(q==1)   
        %---q001の位相でのパケット到着       
        t = statemat(packet+1,1);
        if(~isnan(t))
            A(t,s)=lamdak0+lamdah0+lamdag1;
        else
            loss001 = [loss001 s];
        end
    end

    
   
    if(q==2)
        %---q010の位相でのパケット到着
        t = statemat(packet+1,2);
        if(~isnan(t))
            A(t,s)=lamdak0+lamdah1+lamdag0;
        else
            loss010 = [loss010 s];
        end
    end
    
    
    if(q==3)
        %---q011の位相でのパケット到着
        t = statemat(packet+1,3);
        if(~isnan(t))
            A(t,s)=lamdak0+lamdah1+lamdag1;
        else
            loss011 = [loss011 s];
        end
    end
    
    
    if(q==4)
        %---q100の位相でのパケット到着
        t = statemat(packet+1,4);
        if(~isnan(t))
            A(t,s)=lamdak1+lamdah0+lamdag0;
        else
            loss100 = [loss100 s];
        end
    end
    
    

    if(q==5)
        %---q101の位相でのパケット到着
        t = statemat(packet+1,5);
        if(~isnan(t))
            A(t,s)=lamdak1+lamdah0+lamdag1;
        else
            loss101 = [loss101 s];
        end
    end
    
    

    if(q==6)
        %---q110の位相でのパケット到着
        t = statemat(packet+1,6);
        if(~isnan(t))
            A(t,s)=lamdak1+lamdah1+lamdag0;
        else
            loss110 = [loss110 s];
        end
    end
    

    if(q==7)
        %---q111の位相でのパケット到着
        t = statemat(packet+1,7);
        if(~isnan(t))
            A(t,s)=lamdak1+lamdah1+lamdag1;
        else
            loss111 = [loss111 s];
        end
    end
    
    


    %%----退去

    if(q==0)  %%斜めの遷移を抑制
        t = statemat(packet-1,0);
        if(~isnan(t))        %上のような状態が存在するかチェック．存在しないならt=NaNのはず
            A(t,s)=mu;      %状態がある場合：スタート地点tから行き先sに退去する確率はq000*μ
        end
    end
    
    if(q==1)  %%斜めの遷移を抑制 
        t = statemat(packet-1,1);
        if(~isnan(t))        %上のような状態が存在するかチェック．存在しないならt=NaNのはず
            A(t,s)=mu;      %状態がある場合：スタート地点tから行き先sに退去する確率はq000*μ
        end
    end
    

    if(q==2)  %%斜めの遷移を抑制
        t = statemat(packet-1,2);
        if(~isnan(t))
            A(t,s)=mu;   
        end
    end

    if(q==3)  %%斜めの遷移を抑制 
        t = statemat(packet-1,3);
        if(~isnan(t))
            A(t,s)=mu;   
        end
    end

    if(q==4)  %%斜めの遷移を抑制
        t = statemat(packet-1,4);
        if(~isnan(t))
            A(t,s)=mu;   
        end
    end

    if(q==5)  %%斜めの遷移を抑制 
        t = statemat(packet-1,5);
        if(~isnan(t))
            A(t,s)=mu;
        end
    end

    if(q==6)  %%斜めの遷移を抑制     
        t = statemat(packet-1,6);
        if(~isnan(t))
            A(t,s)=mu;  
        end
    end

    if(q==7)  %%斜めの遷移を抑制 
        t = statemat(packet-1,7);
        if(~isnan(t))
            A(t,s)=mu; 
        end
    end

    %%位相の遷移
    switch q
        case 0
            t = statemat(packet,1);  %%%%%%%+   %%三つの呼のうちどれか1つの位相が変わるパターン
            if(~isnan(t))        %上のような状態が存在するかチェック．存在しないならt=NaNのはず
                A(t,s)=r0g;          %状態がある場合：位相の変化する確率はr0g
            end
            t = statemat(packet,2);
            if(~isnan(t))
                A(t,s)=r0h;
            end
            t = statemat(packet,4);
            if(~isnan(t))
                A(t,s)=r0k; 
            end
        case 1
            t = statemat(packet,0);
            if(~isnan(t))
                A(t,s)=r1g; 
            end
           t = statemat(packet,3);
            if(~isnan(t))
                A(t,s)=r0h; 
            end
            t = statemat(packet,5);
            if(~isnan(t))
                A(t,s)=r0k; 
            end
        case 2
            t = statemat(packet,0);
            if(~isnan(t))
                A(t,s)=r1h; 
            end
            t = statemat(packet,3);
            if(~isnan(t))
                A(t,s)=r0g; 
            end
            t = statemat(packet,6);
            if(~isnan(t))
                A(t,s)=r0k; 
            end
        case 3
            t = statemat(packet,1);
            if(~isnan(t))
                A(t,s)=r1h; 
            end
            t = statemat(packet,2);
            if(~isnan(t))
                A(t,s)=r1g; 
            end
            t = statemat(packet,7);
            if(~isnan(t))
                A(t,s)=r0k; 
            end
        case 4
            t = statemat(packet,0);
            if(~isnan(t))
                A(t,s)=r1k; 
            end
            t = statemat(packet,5);
            if(~isnan(t))
                A(t,s)=r0g; 
            end
            t = statemat(packet,6);
            if(~isnan(t))
                A(t,s)=r0h; 
            end
        case 5
            t = statemat(packet,1);
            if(~isnan(t))
                A(t,s)=r1k; 
            end
            t = statemat(packet,4);
            if(~isnan(t))
                A(t,s)=r0g; 
            end
            t = statemat(packet,7);
            if(~isnan(t))
                A(t,s)=r1h; 
            end
        case 6
            t = statemat(packet,2);
            if(~isnan(t))
                A(t,s)=r1k; 
            end
            t = statemat(packet,4);
            if(~isnan(t))
                A(t,s)=r1h; 
            end
            t = statemat(packet,7);
            if(~isnan(t))
                A(t,s)=r0g; 
            end
        case 7
            t = statemat(packet,3);
            if(~isnan(t))
                A(t,s)=r1k; 
            end
            t = statemat(packet,5);
            if(~isnan(t))
                A(t,s)=r1h; 
            end
            t = statemat(packet,6);
            if(~isnan(t))
                A(t,s)=r1g; 
            end
        otherwise
            fprintf("error\n");
    end

end

%% 各状態の定常確率の計算
P = A+diag(1-sum(A,1));                     %列和が1になるように対角要素をつくる．
[U V] = eigs(P,1);                          %最大固有値に対応する固有ベクトルを計算
Prob = U/sum(U);                             %確率なので，総和が1になるように定数倍



%% パケットロス率の計算
r000 = sum(Prob(loss000));
r001 = sum(Prob(loss001));
r010 = sum(Prob(loss010));
r011 = sum(Prob(loss011));
r100 = sum(Prob(loss100));
r101 = sum(Prob(loss101));
r110 = sum(Prob(loss110));
r111 = sum(Prob(loss111));

totalr = ((r000+r001+r010+r011+r100+r101+r110+r111)/8);
end

