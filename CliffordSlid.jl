using Statistics
using QuantumClifford
using Random
using StableRNGs; 
rng = StableRNG(777);

function ZXZstate(N_sites);
    p0 = zeros(N_sites);p0 = convert.(Bool,p0);
    x1 = copy(p0);x1[1]=1;z1 = copy(p0);z1[2] = 1;z1[end] = 1;
    Pind = PauliOperator(0x0,x1,z1);
    for k in 2:N_sites-1
        x1 = copy(p0);x1[k] = 1;
        z1 = copy(p0);z1[k-1] = 1;z1[k+1] = 1;
        zxz = PauliOperator(0x0,x1,z1);
        Pind = [Pind;zxz]
    end
    x1 = copy(p0);x1[end]=1;z1 = copy(p0);z1[1] = 1;z1[N_sites-1] = 1;
    Pind =[Pind;PauliOperator(0x0,x1,z1)];
    ψ_zxz = Stabilizer(Pind);
    return ψ_zxz
end

function Weight_RC_slide(k)
    w = 1/k/(2.0^k+1)+(k-1)/k/(2.0^k+1)^2;
    return 1/w
end

N_sites = 5;
ψ_ghz = ghz(N_sites);ψ_zxz = ZXZstate(N_sites);
Sample = 1e6;Zexp = [];n0=3;
Omk = [];Ovk = [];Wei_rc = [];k_sub = collect(5:15);Oave = [];
for ks in 1:length(k_sub)
    k_regim = k_sub[ks];
    N_sites = n0 * k_regim;
    ψ_ghz = ghz(N_sites);ψ_zxz = ZXZstate(N_sites);
    state = copy(ψ_zxz);Oave = [];
    n1 = rand(1:N_sites-k_regim+1);n2 = n1+k_regim-1;
    x1 = zeros(N_sites);x1 = convert.(Bool,x1);
    #z1 = copy(x1);z1[n1:n2] .= Bool[1];
    z1 = copy(x1);z1[n1:n1+1] .= 1;z1[n2-1:n2] .= 1;x1[n1+1:n2-1] .= Bool[1]
    Zmeasure = PauliOperator(0x0,x1,z1);
    push!(Zexp,expect(Zmeasure,state));
    Nk = [collect(1:N_sites);collect(1:k_regim)];
    for saa in 1:Sample
        Ac = [];Ai = [];
        s0 = copy(state);
        zk1 = rand(1:k_regim);
        for k1 in 1:Int64(floor(N_sites/k_regim));
            a=random_clifford(rng,k_regim);ai = inv(a);
            push!(Ac,copy(a));push!(Ai,copy(ai));
            zk2 = zk1 + (k1-1)*k_regim;
            apply!(s0,a,Nk[zk2:zk2+k_regim-1]);
        end
        bmix = MixedDestabilizer(s0)
        for k in 1:N_sites
            projectZrand!(bmix,k)
        end
        bmix = stabilizerview(bmix);
        c0 = copy(bmix);
        
        for k1 in 1:Int64(floor(N_sites/k_regim));
            zk2 = zk1 + (k1-1)*k_regim;
            apply!(c0,Ai[k1],Nk[zk2:zk2+k_regim-1]);
        end
        
        Omeas = expect(Zmeasure,c0);
        push!(Oave,Omeas);
    end
    Wrc = Weight_RC_slide(Float64(k_regim));push!(Wei_rc,Wrc);
    Oave .*= Wrc;
    Oest = mean(Oave);
    Ovar = mean(Oave .^2)-Oest.^2;
    push!(Omk,Oest);push!(Ovk,Ovar);
end
