function [C] = FindC(Bq,q,Num)
C = sym(zeros(Num));
for i = 1:Num
    for j = 1:Num
        CIJK = sym([0]);
        for k = 1:Num
            q_dot_temp = q(k)+"_d";
            %fprintf("i = " + i + " j = "+j+" k = "+k+"\n")
            Bij_qk = (diff(Bq(i,j),str2sym(q(k))));
            Bik_qj = (diff(Bq(i,k),str2sym(q(j))));
            Bjk_qi = (diff(Bq(j,k),str2sym(q(i))));
            TempCIJK = 1/2*(Bij_qk + Bik_qj - Bjk_qi);
            
            
            %disp(TempCIJK);
            
            CIJK = CIJK + (TempCIJK*str2sym(q_dot_temp));
            
        end
        C(i,j) = CIJK;
    end
end
C = simplify(C);

end

