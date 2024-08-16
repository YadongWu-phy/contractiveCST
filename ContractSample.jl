using Statistics
using CSV
using DataFrames
using QuantumClifford
using Random
using StableRNGs; 
rng = StableRNG(77);

function ZXZstate(N_sites);
    # ZXZ state
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

function Weight_ZZ(n)
    w = ((5/9)^n-(1/9)^n)/2+(1+(-1/3)^n)/2/3^n;
    return 1/w
end

N_sites = 20;
ψ_ghz = ghz(N_sites);ψ_zxz = ZXZstate(N_sites);
Sample = 1e5;Zexp = [];
Omk = [];Ovk = [];Wei_zz = [];k_sub = collect(5:15);Oave = [];
for ks in 1:length(k_sub)
    k_regim = k_sub[ks];
    ct2 = Contract2();ct2i = inv(ct2);
    state = copy(ψ_zxz);Oave = [];
    n1 = 1;n2 = k_regim;
    x1 = zeros(N_sites);x1 = convert.(Bool,x1);
    #z1 = copy(x1);z1[n1:n2] .= Bool[1];
    #Z₁...Zₖ string operator for GHZ state.
    z1 = copy(x1);z1[n1:n1+1] .= 1;z1[n2-1:n2] .= 1;x1[n1+1:n2-1] .= Bool[1]
    #Z₁Y₂X₃... string operator for ZXZ state.
    Zmeasure = PauliOperator(0x0,x1,z1);
    push!(Zexp,expect(Zmeasure,state));
    for saa in 1:Sample
        a1=[random_clifford1(rng,j) for j in 1:k_regim];ts1 = [CliffordOperator(a, k_regim) for a in a1];
        ti1 = [inv(a) for a in ts1];
        a2=[random_clifford1(rng,j) for j in 1:k_regim];ts2 = [CliffordOperator(a, k_regim) for a in a2];
        ti2 = [inv(a) for a in ts2];
        s0 = copy(state);
        for k1 in 1:k_regim
            apply!(s0,ts1[k1],collect(n1:n2));
        end
        for i1 in n1:n2
            for i2 in i1+1:n2
                apply!(s0,ct2,[i1;i2]);
            end
        end
        for k1 in 1:k_regim
            apply!(s0,ts2[k1],collect(n1:n2));
        end
        bmix = MixedDestabilizer(s0) # |000⟩+|111⟩
        for k in n1:n2
            projectZrand!(bmix,k)
        end
        bmix = stabilizerview(bmix);
        c0 = copy(bmix);
        for k1 in 1:k_regim
            apply!(c0,ti2[k1],collect(n1:n2));
        end
        for i1 in n1:n2
            for i2 in i1+1:n2
                apply!(c0,ct2i,[i1;i2]);
            end
        end
        for k1 in 1:k_regim
            apply!(c0,ti1[k1],collect(n1:n2));
        end
        
        Omeas = expect(Zmeasure,c0);
        push!(Oave,Omeas);
    end
    Wzz = Weight_ZZ(k_regim);push!(Wei_zz,Wzz);# shadow norm for contractive unitary
    Oave .*= Wzz;
    Oest = mean(Oave);
    Oacc = (-1)^k_regim;
    #Oacc = ((-1)^k_regim+1)/2; # GHZ state with ZZZZZ string operator 
    Ovar = mean(Oave .^2)- Oest^2;
    push!(Omk,Oest);push!(Ovk,Ovar);
end
nfile = join(["N",string(N_sites),"Contract_ZXZ.CSV"])
name1 = ["N_sites","k_sub","Omk","Ovk","W_ct","Msample"];
valu1 = [N_sites,k_sub,Omk,Ovk,Wzz,Sample];
df1 = DataFrame(name = name1;value = valu1);
CSV.write(nfile,df1)
