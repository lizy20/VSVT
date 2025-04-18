function [ke,D]=kemefe(nu,le)
E0 =1;
D=E0/(1-nu^2)*[1 nu 0;
    nu 1 0;
    0 0 (1-nu)/2];
ke=zeros(8,8);
NODEXY=[0,0;le,0;le,le;0,le]; % 瀹炲弬鍧愭爣鍗曞厓鍧愭爣 鍙煡锛歔杈归暱澶у皬]

sq=sqrt(3)/3;
l=[-sq -sq; sq -sq; sq sq; -sq sq];
w=[1;1;1;1];

for q=1:4
    gausspoint=l(q,:);
    xi=gausspoint(1);
    eta=gausspoint(2);

    % Shape functions and derivatives
    N=[1/4*(1-xi)*(1-eta) 1/4*(1+xi)*(1-eta) 1/4*(1+xi)*(1+eta) 1/4*(1-xi)*(1+eta)]; % 褰㈢姸鍑芥暟
    dN=1/4*[-(1-eta) 1-eta  1+eta -(1+eta); 
            -(1-xi) -(1+xi) 1+xi  1-xi]; % 褰㈠嚱鏁板鏁?

    % Jacobian matrix, inverse of Jacobian
    J=dN*NODEXY; % 闆呭彲姣旂煩闃? = 璁＄畻鏀惧ぇ鍊嶆暟 锛岃嚜鐒跺潗鏍囩郴涓嬬殑澶у皬锛氬眬閮ㄥ潗鏍囩郴涓嬬殑澶у皬
    dNxy=inv(J)*dN; % dxdy = det[J]*d尉d畏 灏嗗眬閮ㄥ潗鏍囩郴涓嬬殑绉垎鎹㈢畻鍒拌嚜鐒跺潗鏍囩郴涓?

    % B matrix
    B=[dNxy(1,1) 0 dNxy(1,2) 0 dNxy(1,3) 0 dNxy(1,4) 0;
        0 dNxy(2,1) 0 dNxy(2,2) 0 dNxy(2,3) 0 dNxy(2,4);
        dNxy(2,1) dNxy(1,1) dNxy(2,2) dNxy(1,2) dNxy(2,3) dNxy(1,3) dNxy(2,4) dNxy(1,4)];
    % stiffness matri
    ke=ke+B'*D*B*w(q)*det(J);
    

end
end