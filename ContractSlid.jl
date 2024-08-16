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

function Contract2()
    p0 = zeros(2);p0 = convert.(Bool,p0);
    z0 = ones(2); z0 = convert.(Bool,z0);
    x1 = copy(p0);x1[1]=1;
    Pind = PauliOperator(0x2,x1,z0);
    x1 = copy(p0);x1[2] = 1;
    zxz = PauliOperator(0x2,x1,z0);
    Pind = [Pind;zxz]
    for k in 1:2
        z1 = copy(p0);z1[k] = 1;
        zxz = PauliOperator(0x0,p0,z1);
        Pind = [Pind;zxz]
    end
    Cli_zz = CliffordOperator(Pind);
    return Cli_zz;
end
    
function Weight_ZZ_slide(k)
    w = 1/64*27^(-k)*(24+4*(-3)^k-4*(-3)^(2*k)-6*(-1)^k
        -4*5^k+4*3^k-4*(-1)^k*3^(1+k)+5^k*3^(1+k)-4*3^(2+k)
        +4*3^(1+2*k)+5^(1+k)*(-4+7*3^k)+16*(1+(-3)^k-2*3^k+5^k+9^k)*k);
    return k/w
end

N_sites = 10;
ψ_ghz = ghz(N_sites);ψ_zxz = ZXZstate(N_sites);
Sample = 1e6;Zexp = [];n0 = 2;
Omk = [];Ovk = [];Wei_zz = [];k_sub = collect(5:15);Oave = [];
for ks in 1:length(k_sub)
    k_regim = k_sub[ks];
    N_sites = n0*k_regim;
    ψ_ghz = ghz(N_sites);ψ_zxz = ZXZstate(N_sites);
    ct2 = Contract2();ct2i = inv(Cli_zz);
    state = copy(ψ_zxz);Oave = [];
    n1 = rand(1:N_sites-k_regim+1);n2 = n1+k_regim-1;
    x1 = zeros(N_sites);x1 = convert.(Bool,x1);
    #z1 = copy(x1);z1[n1:n2] .= Bool[1];
    z1 = copy(x1);z1[n1:n1+1] .= 1;z1[n2-1:n2] .= 1;x1[n1+1:n2-1] .= Bool[1]
    Zmeasure = PauliOperator(0x0,x1,z1);
    push!(Zexp,expect(Zmeasure,state));
    Nk = [collect(1:N_sites);collect(1:k_regim)];
    for saa in 1:Sample
        a1=[random_clifford1(rng,j) for j in 1:N_sites];ts1 = [CliffordOperator(a, N_sites) for a in a1];
        ti1 = [inv(a) for a in ts1];
        a2=[random_clifford1(rng,j) for j in 1:N_sites];ts2 = [CliffordOperator(a, N_sites) for a in a2];
        ti2 = [inv(a) for a in ts2];
        s0 = copy(state);
        for k1 in 1:N_sites
            apply!(s0,ts1[k1],collect(1:N_sites));
        end
        zk1 = rand(1:k_regim);
        for k1 in 1:Int64(floor(N_sites/k_regim));
            zk2 = zk1 + (k1-1)*k_regim;
            nreg = Nk[zk2:zk2+k_regim-1];
            for i1 in 1:k_regim
                for i2 in i1+1:k_regim
                    apply!(s0,ct2,[nreg[i1];nreg[i2]]);
                end
            end
        end
        for k1 in 1:N_sites
            apply!(s0,ts2[k1],collect(1:N_sites));
        end
        bmix = MixedDestabilizer(s0) # |000⟩+|111⟩
        for k in 1:N_sites
            projectZrand!(bmix,k)
        end
        bmix = stabilizerview(bmix);
        c0 = copy(bmix);
        for k1 in 1:N_sites
            apply!(c0,ti2[k1],collect(1:N_sites));
        end
        for k1 in 1:Int64(floor(N_sites/k_regim));
            zk2 = zk1 + (k1-1)*k_regim;
            nreg = Nk[zk2:zk2+k_regim-1];
            for i1 in 1:k_regim
                for i2 in i1+1:k_regim
                    apply!(c0,ct2i,[nreg[i1];nreg[i2]]);
                end
            end
        end
        for k1 in 1:N_sites
            apply!(c0,ti1[k1],collect(1:N_sites));
        end
        
        Omeas = expect(Zmeasure,c0);
        push!(Oave,Omeas);
    end
    Wzz = Weight_ZZ_slide(Float64(k_regim));push!(Wei_zz,Wzz);
    Oave .*= Wzz;
    Oest = mean(Oave);
    Ovar = mean(Oave .^2)-Oest.^2;
    push!(Omk,Oest);push!(Ovk,Ovar);
end
