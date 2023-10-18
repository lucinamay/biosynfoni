leaf = [
    {"name": "Isobutyric acid", "smarts": "[CH3][CH1]([CH3])[CH0]([OH])=[O]"},
    {
        "name": "k-Arg",
        "smarts": "[NH]=[C]([NH2])[NH][CH2][CH2][CH2][CH]([NH2])[C](=[O])[C](=[O])[OH]",
    },
    {
        "name": "Ac-Trp",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][c]1[cH][nH][c]2[cH][cH][cH][cH][c]12)[C](=[O])[OH]",
    },
    {
        "name": "iC13:1(3)",
        "smarts": "[CH3][CH]([CH3])[CH2][CH2][CH2][CH2][CH2][CH2][CH]=[CH][CH2][C](=[O])[OH]",
    },
    {
        "name": "iC16:0-OH(3)",
        "smarts": "[CH3][CH]([CH3])[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH]([OH])[CH2][C](=[O])[OH]",
    },
    {"name": "bAbu", "smarts": "[CH3][CH]([NH2])[CH2][C](=[O])[OH]"},
    {
        "name": "3Me-Phe",
        "smarts": "[CH3][c]1[cH][cH][cH][c]([CH2][CH]([NH2])[C](=[O])[OH])[cH]1",
    },
    {
        "name": "C6:0-Me(5.5)-oxo(2)",
        "smarts": "[CH3][C]([CH3])([CH3])[CH2][CH2][C](=[O])[C](=[O])[OH]",
    },
    {"name": "NMe-Abu", "smarts": "[CH3][CH2][CH]([NH][CH3])[C](=[O])[OH]"},
    {
        "name": "NAc-4OH-Pro",
        "smarts": "[CH3][C](=[O])[N]1[CH2][CH]([OH])[CH2][CH]1[C](=[O])[OH]",
    },
    {
        "name": "NMe-t-Leu",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[C]([CH3])([CH3])[CH3]",
    },
    {"name": "OH-4Abu", "smarts": "[NH2][CH2][CH]([OH])[CH2][C](=[O])[OH]"},
    {
        "name": "NAc-Dbu",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[CH]([CH3])[NH2]",
    },
    {"name": "Spd", "smarts": "[NH2][CH2][CH2][CH2][CH2][NH][CH2][CH2][CH2][NH2]"},
    {
        "name": "Me-AOA",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH]([NH2])[C](=[O])[OH]",
    },
    {
        "name": "Br-OH-Trp",
        "smarts": "[NH2][CH]([CH2][c]1[c]([Br])[nH][c]2[cH][cH][c]([OH])[cH][c]12)[C](=[O])[OH]",
    },
    {
        "name": "3OH-5Me-Pro",
        "smarts": "[CH3][CH]1[CH2][CH]([OH])[CH]([C](=[O])[OH])[NH]1",
    },
    {
        "name": "C12:0",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][C](=[O])[OH]",
    },
    {
        "name": "NMe-3OH-5Me-Pro",
        "smarts": "[CH3][CH]1[CH2][CH]([OH])[CH]([C](=[O])[OH])[N]1[CH3]",
    },
    {
        "name": "NFo-OMe-Trp",
        "smarts": "[CH3][O][n]1[cH][c]([CH2][CH]([NH][CH]=[O])[C](=[O])[OH])[c]2[cH][cH][cH][cH][c]21",
    },
    {
        "name": "NAc-Hpr",
        "smarts": "[CH3][C](=[O])[N]1[CH2][CH2][CH2][CH2][CH]1[C](=[O])[OH]",
    },
    {
        "name": "N1-COOH-bhTrp",
        "smarts": "[NH2][CH]([CH2][CH2][CH2][c]1[cH][n]([C](=[O])[OH])[c]2[cH][cH][cH][cH][c]12)[C](=[O])[OH]",
    },
    {"name": "b4Cl-Thr", "smarts": "[NH2][C]([OH])([CH2][Cl])[CH2][C](=[O])[OH]"},
    {
        "name": "NFo-bOH-Cl-Tyr",
        "smarts": "[O]=[CH][NH][CH]([C](=[O])[OH])[CH]([OH])[c]1[cH][cH][c]([OH])[c]([Cl])[cH]1",
    },
    {
        "name": "NAc-Hil",
        "smarts": "[CH3][CH2][CH]([CH3])[CH2][CH]([NH][C]([CH3])=[O])[C](=[O])[OH]",
    },
    {
        "name": "NMe-Har",
        "smarts": "[CH3][NH][CH]([CH2][CH2][CH2][CH2][N]=[C]([NH2])[NH2])[C](=[O])[OH]",
    },
    {"name": "NMe-Thr", "smarts": "[CH3][NH][CH]([C](=[O])[OH])[CH]([CH3])[OH]"},
    {
        "name": "bbMe2-O-Met",
        "smarts": "[CH3][S](=[O])[CH2][C]([CH3])([CH3])[CH]([NH2])[C](=[O])[OH]",
    },
    {"name": "bAla", "smarts": "[NH2][CH2][CH2][C](=[O])[OH]"},
    {
        "name": "bHph",
        "smarts": "[NH2][CH]([CH2][C](=[O])[OH])[CH2][c]1[cH][cH][cH][cH][cH]1",
    },
    {"name": "C6:0-Ep(2)", "smarts": "[CH3][CH2][CH2][CH]1[O][CH]1[C](=[O])[OH]"},
    {
        "name": "bOH-Tyr",
        "smarts": "[NH2][CH]([C](=[O])[OH])[CH]([OH])[c]1[cH][cH][c]([OH])[cH][cH]1",
    },
    {
        "name": "NMe-C10:0-OH(9)-NH2(2)",
        "smarts": "[CH3][NH][CH]([CH2][CH2][CH2][CH2][CH2][CH2][CH]([CH3])[OH])[C](=[O])[OH]",
    },
    {
        "name": "NFo-NMe-Lan",
        "smarts": "[O]=[CH][NH][CH]([CH2][S][CH2][CH]([NH][OH])[C](=[O])[OH])[C](=[O])[OH]",
    },
    {
        "name": "diOH-Arg",
        "smarts": "[NH2][C]([NH2])=[N][CH2][CH]([OH])[CH]([OH])[CH]([NH2])[C](=[O])[OH]",
    },
    {
        "name": "C8:0-Me(4)-OH(3)",
        "smarts": "[CH3][CH2][CH2][CH2][CH]([CH3])[CH]([OH])[CH2][C](=[O])[OH]",
    },
    {
        "name": "NFo-3Me-Hty",
        "smarts": "[CH3][CH]([CH2][c]1[cH][cH][c]([OH])[cH][cH]1)[CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "NAc-Br-Trp",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][c]1[cH][nH][c]2[cH][cH][c]([Br])[cH][c]12)[C](=[O])[OH]",
    },
    {
        "name": "C11:0-Me(2)-OH(3)",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH]([OH])[CH]([CH3])[C](=[O])[OH]",
    },
    {
        "name": "NFo-Ac-Ser",
        "smarts": "[CH3][C](=[O])[O][CH2][CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "NFo-bMe-Phe",
        "smarts": "[CH3][CH]([c]1[cH][cH][cH][cH][cH]1)[CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "MeOx-Ile",
        "smarts": "[CH3][CH2][CH]([CH3])[CH]([NH2])[C]1=[N][CH]([C](=[O])[OH])[CH]([CH3])[O]1",
    },
    {
        "name": "bMe-Gln",
        "smarts": "[CH3][CH]([CH2][C]([NH2])=[O])[CH]([NH2])[C](=[O])[OH]",
    },
    {
        "name": "aC6:0-OH(2.3)",
        "smarts": "[CH3][CH2][C]([CH3])([OH])[CH]([OH])[C](=[O])[OH]",
    },
    {
        "name": "NAc-C10:0-OH(8)-NH2(2)",
        "smarts": "[CH3][CH2][CH]([OH])[CH2][CH2][CH2][CH2][CH2][CH]([NH][C]([CH3])=[O])[C](=[O])[OH]",
    },
    {
        "name": "NAc-Cl-Trp",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][C]1=[N][c]2[cH][c]([Cl])[cH][cH][c]2[CH2]1)[C](=[O])[OH]",
    },
    {
        "name": "NAc-bMe-Br-Phe",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[CH]([CH3])[c]1[cH][cH][c]([Br])[cH][cH]1",
    },
    {
        "name": "NMe-Bmt",
        "smarts": "[CH3][CH]=[CH][CH2][CH]([CH3])[CH]([OH])[CH]([NH][CH3])[C](=[O])[OH]",
    },
    {
        "name": "NFo-t-Leu",
        "smarts": "[CH3][C]([CH3])([CH3])[CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "bC10:0-OH(9)-NH2(2)-oxo(8)",
        "smarts": "[CH3][CH]([OH])[C](=[O])[CH2][CH2][CH2][CH2][CH]([NH2])[CH2][C](=[O])[OH]",
    },
    {
        "name": "2Me-3Me-pGlu",
        "smarts": "[CH3][CH]1[C](=[O])[NH][CH]([C](=[O])[OH])[CH]1[CH3]",
    },
    {
        "name": "NAc-2Me-3Me-pGlu",
        "smarts": "[CH3][C](=[O])[N]1[C](=[O])[CH]([CH3])[CH]([CH3])[CH]1[C](=[O])[OH]",
    },
    {"name": "NAc-Dpr", "smarts": "[CH3][C](=[O])[NH][CH]([CH2][NH2])[C](=[O])[OH]"},
    {
        "name": "NAc-bOMe-Tyr",
        "smarts": "[CH3][O][CH]([c]1[cH][cH][c]([OH])[cH][cH]1)[CH]([NH][C]([CH3])=[O])[C](=[O])[OH]",
    },
    {
        "name": "Ph-Ser",
        "smarts": "[NH2][CH]([C](=[O])[OH])[CH]([OH])[c]1[cH][cH][cH][cH][cH]1",
    },
    {
        "name": "aC9:0",
        "smarts": "[CH3][CH2][CH]([CH3])[CH2][CH2][CH2][CH2][C](=[O])[OH]",
    },
    {"name": "HseL", "smarts": "[NH2][CH]1[CH2][CH2][O][C]1=[O]"},
    {
        "name": "NAc-Ac-Ser",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][O][C]([CH3])=[O])[C](=[O])[OH]",
    },
    {
        "name": "NAc-N2Me-bOH-Asn",
        "smarts": "[CH3][NH][C](=[O])[CH]([OH])[CH]([NH][C]([CH3])=[O])[C](=[O])[OH]",
    },
    {
        "name": "NAc-t-Leu",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[C]([CH3])([CH3])[CH3]",
    },
    {"name": "bOH-Val", "smarts": "[CH3][C]([CH3])([OH])[CH]([NH2])[C](=[O])[OH]"},
    {
        "name": "C10:0-OH(9)-NH2(2)",
        "smarts": "[CH3][CH]([OH])[CH2][CH2][CH2][CH2][CH2][CH2][CH]([NH2])[C](=[O])[OH]",
    },
    {
        "name": "aC9:0-OH(3)",
        "smarts": "[CH3][CH2][CH]([CH3])[CH2][CH2][CH]([OH])[CH2][C](=[O])[OH]",
    },
    {"name": "Hil", "smarts": "[CH3][CH2][CH]([CH3])[CH2][CH]([NH2])[C](=[O])[OH]"},
    {"name": "NFo-Iser", "smarts": "[O]=[CH][NH][CH2][CH]([OH])[C](=[O])[OH]"},
    {
        "name": "NFo-bOH-Val",
        "smarts": "[CH3][C]([CH3])([OH])[CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "bC10:0-NH2(2)-oxo(8)",
        "smarts": "[CH3][CH2][C](=[O])[CH2][CH2][CH2][CH2][CH]([NH2])[CH2][C](=[O])[OH]",
    },
    {
        "name": "NMe-Hph",
        "smarts": "[CH3][NH][CH]([CH2][CH2][c]1[cH][cH][cH][cH][cH]1)[C](=[O])[OH]",
    },
    {"name": "Dov", "smarts": "[CH3][CH]([CH3])[CH]([C](=[O])[OH])[N]([CH3])[CH3]"},
    {"name": "NdMe-Ala", "smarts": "[CH3][CH]([C](=[O])[OH])[N]([CH3])[CH3]"},
    {
        "name": "NMe-bOH-Gln",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[CH]([OH])[CH2][C]([NH2])=[O]",
    },
    {
        "name": "Har",
        "smarts": "[NH2][C]([NH2])=[N][CH2][CH2][CH2][CH2][CH]([NH2])[C](=[O])[OH]",
    },
    {
        "name": "NAc-Cl2-Hpg",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[c]1[cH][c]([Cl])[c]([OH])[c]([Cl])[cH]1",
    },
    {
        "name": "NAc-Asn",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][C]([NH2])=[O])[C](=[O])[OH]",
    },
    {
        "name": "NMe-OH-Orn",
        "smarts": "[CH3][NH][CH]([CH2][CH2][CH2][NH][OH])[C](=[O])[OH]",
    },
    {
        "name": "aC15:0-OH(3)",
        "smarts": "[CH3][CH2][CH]([CH3])[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH]([OH])[CH2][C](=[O])[OH]",
    },
    {
        "name": "OH-His",
        "smarts": "[NH2][CH]([C](=[O])[OH])[CH]([OH])[c]1[cH][n][cH][nH]1",
    },
    {
        "name": "NMe-NAc-Ile",
        "smarts": "[CH3][CH2][CH]([CH3])[CH]([C](=[O])[OH])[N]([CH3])[C]([CH3])=[O]",
    },
    {
        "name": "NFo-bMe-Asn",
        "smarts": "[CH3][CH]([C]([NH2])=[O])[CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "bAhad",
        "smarts": "[NH2][CH]([CH2][C](=[O])[OH])[CH]([OH])[CH2][C](=[O])[OH]",
    },
    {
        "name": "bN1-COOH-bhTrp",
        "smarts": "[NH2][CH]([CH2][CH2][c]1[cH][n]([C](=[O])[OH])[c]2[cH][cH][cH][cH][c]12)[CH2][C](=[O])[OH]",
    },
    {
        "name": "NMe-aIle/NMe-Ile",
        "smarts": "[CH3][CH2][CH]([CH3])[CH]([NH][CH3])[C](=[O])[OH]",
    },
    {"name": "NMe-Gly", "smarts": "[CH3][NH][CH2][C](=[O])[OH]"},
    {
        "name": "NAc-Apv",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][CH2][CH2][c]1[cH][cH][cH][cH][cH]1)[C](=[O])[OH]",
    },
    {
        "name": "OH-dHpg",
        "smarts": "[O]=[C]1[CH]=[CH][C](=[C]([NH][OH])[C](=[O])[OH])[CH]=[CH]1",
    },
    {
        "name": "aC9:0-OH(2)-NH2(3)",
        "smarts": "[CH3][CH2][CH]([CH3])[CH2][CH2][CH]([NH2])[CH]([OH])[C](=[O])[OH]",
    },
    {
        "name": "bHty",
        "smarts": "[NH2][CH]([CH2][C](=[O])[OH])[CH2][c]1[cH][cH][c]([OH])[cH][cH]1",
    },
    {
        "name": "NAc-PO-Asn",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[CH]([O][P](=[O])([OH])[OH])[C]([NH2])=[O]",
    },
    {
        "name": "NAc-OH-Asn",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[CH]([OH])[C]([NH2])=[O]",
    },
    {
        "name": "C11:2(t2.t8)-Me(2.6.8)-OH(5.7)",
        "smarts": "[CH3][CH2][CH]=[C]([CH3])[CH]([OH])[CH]([CH3])[CH]([OH])[CH2][CH]=[C]([CH3])[C](=[O])[OH]",
    },
    {
        "name": "NFo-bMe-Asp",
        "smarts": "[CH3][CH]([C](=[O])[OH])[CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {"name": "Hiv", "smarts": "[CH3][CH]([CH3])[CH]([OH])[C](=[O])[OH]"},
    {
        "name": "NMe-Ac-OH-Orn",
        "smarts": "[CH3][NH][CH]([CH2][CH2][CH2][N]([OH])[C]([CH3])=[O])[C](=[O])[OH]",
    },
    {
        "name": "v-Tyr",
        "smarts": "[NH2][CH]([CH]=[CH][C](=[O])[OH])[CH2][c]1[cH][cH][c]([OH])[cH][cH]1",
    },
    {
        "name": "bDaz",
        "smarts": "[NH2][CH]([CH2][CH2][CH]([NH2])[CH]([OH])[CH2][C](=[O])[OH])[CH2][C](=[O])[OH]",
    },
    {
        "name": "Ac-Val",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[CH]([CH3])[CH3]",
    },
    {
        "name": "Trpol",
        "smarts": "[NH2][CH]([CH2][OH])[CH2][c]1[cH][nH][c]2[cH][cH][cH][cH][c]12",
    },
    {
        "name": "aC15:0-NH2(3)",
        "smarts": "[CH3][CH2][CH]([CH3])[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH]([NH2])[CH2][C](=[O])[OH]",
    },
    {
        "name": "NAc-bOH-Tyr",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[CH]([OH])[c]1[cH][cH][c]([OH])[cH][cH]1",
    },
    {
        "name": "NFo-Br-Phe",
        "smarts": "[O]=[CH][NH][CH]([CH2][c]1[cH][cH][c]([Br])[cH][cH]1)[C](=[O])[OH]",
    },
    {
        "name": "C12:3(7.9.11)-Me(6.10)-OH(2.4.5)-NH2(3)-Ph(12)",
        "smarts": "[CH3][C]([CH]=[CH][c]1[cH][cH][cH][cH][cH]1)=[CH][CH]=[CH][CH]([CH3])[CH]([OH])[CH]([OH])[CH]([NH2])[CH]([OH])[C](=[O])[OH]",
    },
    {
        "name": "Ph-Lac",
        "smarts": "[O]=[C]([OH])[CH]([OH])[CH2][c]1[cH][cH][cH][cH][cH]1",
    },
    {"name": "Dpr", "smarts": "[NH2][CH2][CH]([NH2])[C](=[O])[OH]"},
    {
        "name": "NAc-bMe-Asn",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[CH]([CH3])[C]([NH2])=[O]",
    },
    {"name": "NMe-pGlu", "smarts": "[CH3][N]1[C](=[O])[CH2][CH2][CH]1[C](=[O])[OH]"},
    {
        "name": "bAca",
        "smarts": "[NH2][CH]([CH2][C](=[O])[OH])[CH]1[CH2][CH2][C](=[O])[CH]2[O][CH]21",
    },
    {"name": "Serol", "smarts": "[NH2][CH]([CH2][OH])[CH2][OH]"},
    {
        "name": "NFo-N1-COOH-bhTrp",
        "smarts": "[O]=[CH][NH][CH]([CH2][CH2][CH2][c]1[cH][n]([C](=[O])[OH])[c]2[cH][cH][cH][cH][c]12)[C](=[O])[OH]",
    },
    {
        "name": "OMe-bAla-Thz",
        "smarts": "[CH3][O][C](=[O])[CH2][CH]([NH2])[c]1[n][cH][cH][s]1",
    },
    {
        "name": "NFo-C10:0-OH(9)-NH2(2)-oxo(8)",
        "smarts": "[CH3][CH]([OH])[C](=[O])[CH2][CH2][CH2][CH2][CH2][CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "C16:0-OH(3.4)",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH]([OH])[CH]([OH])[CH2][C](=[O])[OH]",
    },
    {
        "name": "NFo-2Me-3Me-pGlu",
        "smarts": "[CH3][CH]1[C](=[O])[N]([CH]=[O])[CH]([C](=[O])[OH])[CH]1[CH3]",
    },
    {"name": "NFo-pGlu", "smarts": "[O]=[CH][N]1[C](=[O])[CH2][CH2][CH]1[C](=[O])[OH]"},
    {
        "name": "ADMAdda",
        "smarts": "[CH3][C](=[O])[O][CH]([CH2][c]1[cH][cH][cH][cH][cH]1)[CH]([CH3])[CH]=[C]([CH3])[CH]=[CH][CH]([NH2])[CH]([CH3])[C](=[O])[OH]",
    },
    {"name": "C4:0-OH(3)", "smarts": "[CH3][CH]([OH])[CH2][C](=[O])[OH]"},
    {"name": "bN2Me-Asn", "smarts": "[CH3][NH][C](=[O])[CH]([NH2])[CH2][C](=[O])[OH]"},
    {
        "name": "C10:0",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][C](=[O])[OH]",
    },
    {"name": "iC5:0-OH(2.4)", "smarts": "[CH3][CH]([CH2][OH])[CH]([OH])[C](=[O])[OH]"},
    {
        "name": "NMe-Br-Phe",
        "smarts": "[CH3][NH][CH]([CH2][c]1[cH][cH][c]([Br])[cH][cH]1)[C](=[O])[OH]",
    },
    {
        "name": "NFo-Dhpg",
        "smarts": "[O]=[CH][NH][CH]([C](=[O])[OH])[c]1[cH][c]([OH])[cH][c]([OH])[cH]1",
    },
    {
        "name": "NFo-4Cl-Thr",
        "smarts": "[O]=[CH][NH][CH]([C](=[O])[OH])[CH]([OH])[CH2][Cl]",
    },
    {
        "name": "bFo-OH-Orn",
        "smarts": "[NH2][CH]([CH2][CH2][N]([OH])[CH]=[O])[CH2][C](=[O])[OH]",
    },
    {
        "name": "NFo-Kyn",
        "smarts": "[NH2][c]1[cH][cH][cH][cH][c]1[C](=[O])[CH2][CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {"name": "norCMA", "smarts": "[CH3][CH]1[CH2][C]1([NH2])[C](=[O])[OH]"},
    {"name": "C4:0", "smarts": "[CH3][CH2][CH2][C](=[O])[OH]"},
    {
        "name": "Dil",
        "smarts": "[CH3][CH2][CH]([CH3])[CH]([NH][CH3])[CH]([CH2][C](=[O])[OH])[O][CH3]",
    },
    {
        "name": "OH-bLys",
        "smarts": "[NH2][CH2][CH2][CH]([OH])[CH]([NH2])[CH2][C](=[O])[OH]",
    },
    {
        "name": "Cl3-NMe-Leu",
        "smarts": "[CH3][NH][CH]([CH2][CH]([CH3])[C]([Cl])([Cl])[Cl])[C](=[O])[OH]",
    },
    {
        "name": "Choi",
        "smarts": "[O]=[C]([OH])[CH]1[CH2][CH]2[CH2][CH2][CH]([OH])[CH2][CH]2[NH]1",
    },
    {
        "name": "iC15:0-NH2(3)",
        "smarts": "[CH3][CH]([CH3])[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH]([NH2])[CH2][C](=[O])[OH]",
    },
    {"name": "Asn", "smarts": "[NH2][C](=[O])[CH2][CH]([NH2])[C](=[O])[OH]"},
    {
        "name": "NFo-CysA",
        "smarts": "[O]=[CH][NH][CH]([CH2][S](=[O])(=[O])[OH])[C](=[O])[OH]",
    },
    {
        "name": "NFo-Hph",
        "smarts": "[O]=[CH][NH][CH]([CH2][CH2][c]1[cH][cH][cH][cH][cH]1)[C](=[O])[OH]",
    },
    {
        "name": "1Me-Trp",
        "smarts": "[CH3][n]1[cH][c]([CH2][CH]([NH2])[C](=[O])[OH])[c]2[cH][cH][cH][cH][c]21",
    },
    {
        "name": "NAc-F-ph-Gly",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[c]1[cH][cH][c]([F])[cH][cH]1",
    },
    {
        "name": "NMe-OH-His",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[CH]([OH])[c]1[cH][n][cH][nH]1",
    },
    {"name": "bEt-Nva", "smarts": "[CH3][CH2][C]([NH2])([CH2][CH3])[CH2][C](=[O])[OH]"},
    {
        "name": "Tyr",
        "smarts": "[NH2][CH]([CH2][c]1[cH][cH][c]([OH])[cH][cH]1)[C](=[O])[OH]",
    },
    {
        "name": "NMe-C10:0-OH(9)-NH2(2)-oxo(8)",
        "smarts": "[CH3][NH][CH]([CH2][CH2][CH2][CH2][CH2][C](=[O])[CH]([CH3])[OH])[C](=[O])[OH]",
    },
    {"name": "Hap", "smarts": "[CH3][CH]([C](=[O])[OH])[C](=[O])[CH2][OH]"},
    {
        "name": "NMe-hv-Val",
        "smarts": "[CH3][NH][CH]([CH]=[C]([CH3])[C](=[O])[OH])[CH]([CH3])[CH3]",
    },
    {"name": "Abu", "smarts": "[CH3][CH2][CH]([NH2])[C](=[O])[OH]"},
    {
        "name": "iC9:2(2.t4)",
        "smarts": "[CH3][CH]([CH3])[CH2][CH]=[CH][CH]=[CH][C](=[O])[OH]",
    },
    {
        "name": "NFo-OH-Asp",
        "smarts": "[O]=[CH][NH][CH]([C](=[O])[OH])[CH]([OH])[C](=[O])[OH]",
    },
    {
        "name": "O2-Met",
        "smarts": "[CH3][S](=[O])(=[O])[CH2][CH2][CH]([NH2])[C](=[O])[OH]",
    },
    {
        "name": "NAc-3Me-Glu",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[CH]([CH3])[CH2][C](=[O])[OH]",
    },
    {
        "name": "NMe-Hty",
        "smarts": "[CH3][NH][CH]([CH2][CH2][c]1[cH][cH][c]([OH])[cH][cH]1)[C](=[O])[OH]",
    },
    {
        "name": "NFo-C10:0-NH2(2)-oxo(8)",
        "smarts": "[CH3][CH2][C](=[O])[CH2][CH2][CH2][CH2][CH2][CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "NFo-C10:0-OH(9)-NH2(2)",
        "smarts": "[CH3][CH]([OH])[CH2][CH2][CH2][CH2][CH2][CH2][CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "3d-NMe-Bmt",
        "smarts": "[CH3][CH]=[CH][CH2][CH]([CH3])[CH2][CH]([NH][CH3])[C](=[O])[OH]",
    },
    {
        "name": "C16:1(9)",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH]=[CH][CH2][CH2][CH2][CH2][CH2][CH2][CH2][C](=[O])[OH]",
    },
    {
        "name": "Br-Tyr",
        "smarts": "[NH2][CH]([CH2][c]1[cH][cH][c]([OH])[c]([Br])[cH]1)[C](=[O])[OH]",
    },
    {"name": "NFo-Nva", "smarts": "[CH3][CH2][CH2][CH]([NH][CH]=[O])[C](=[O])[OH]"},
    {
        "name": "NMe-Me2A-Phe",
        "smarts": "[CH3][NH][CH]([CH2][c]1[cH][cH][c]([N]([CH3])[CH3])[cH][cH]1)[C](=[O])[OH]",
    },
    {
        "name": "NFo-bbMe2-O-Met",
        "smarts": "[CH3][S](=[O])[CH2][C]([CH3])([CH3])[CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "iC14:0-OH(3)",
        "smarts": "[CH3][CH]([CH3])[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH]([OH])[CH2][C](=[O])[OH]",
    },
    {
        "name": "C12:0-OH(3)",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH]([OH])[CH2][C](=[O])[OH]",
    },
    {"name": "aFo-Gly", "smarts": "[NH2][CH]([CH]=[O])[C](=[O])[OH]"},
    {"name": "MdCP", "smarts": "[CH3][n]1[c]([C](=[O])[OH])[cH][c]([Cl])[c]1[Cl]"},
    {
        "name": "NMe-OH-Ile",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[CH]([CH3])[CH2][CH2][OH]",
    },
    {
        "name": "C16:1(9)-OH(3)",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH]=[CH][CH2][CH2][CH2][CH2][CH2][CH]([OH])[CH2][C](=[O])[OH]",
    },
    {
        "name": "C14:0-OH(3)",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH]([OH])[CH2][C](=[O])[OH]",
    },
    {
        "name": "NAc-bMe-Gln",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[CH]([CH3])[CH2][C]([NH2])=[O]",
    },
    {
        "name": "5OH-Cap",
        "smarts": "[NH]=[C]1[NH][CH]([OH])[CH2][CH]([CH]([NH2])[C](=[O])[OH])[NH]1",
    },
    {
        "name": "NAc-Nva",
        "smarts": "[CH3][CH2][CH2][CH]([NH][C]([CH3])=[O])[C](=[O])[OH]",
    },
    {
        "name": "NAc-OH-Asp",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[CH]([OH])[C](=[O])[OH]",
    },
    {
        "name": "NFo-Aad",
        "smarts": "[O]=[CH][NH][CH]([CH2][CH2][CH2][C](=[O])[OH])[C](=[O])[OH]",
    },
    {
        "name": "NMe-Lan",
        "smarts": "[NH2][CH]([CH2][S][CH2][CH]([NH][OH])[C](=[O])[OH])[C](=[O])[OH]",
    },
    {
        "name": "NFo-Har",
        "smarts": "[NH2][C]([NH2])=[N][CH2][CH2][CH2][CH2][CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "NAc-4Me-Pro",
        "smarts": "[CH3][C](=[O])[N]1[CH2][CH]([CH3])[CH2][CH]1[C](=[O])[OH]",
    },
    {"name": "bCl-Ile", "smarts": "[CH3][CH]([Cl])[C]([CH3])([NH2])[CH2][C](=[O])[OH]"},
    {
        "name": "NMe-5Me-Pro",
        "smarts": "[CH3][CH]1[CH2][CH2][CH]([C](=[O])[OH])[N]1[CH3]",
    },
    {
        "name": "NAc-Dhpg",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[c]1[cH][c]([OH])[cH][c]([OH])[cH]1",
    },
    {
        "name": "NFo-4oxo-5Me-Pro",
        "smarts": "[CH3][CH]1[C](=[O])[CH2][CH]([C](=[O])[OH])[N]1[CH]=[O]",
    },
    {
        "name": "Dap",
        "smarts": "[CH3][O][CH]([CH]1[CH2][CH2][CH2][NH]1)[CH]([CH3])[C](=[O])[OH]",
    },
    {"name": "Dab", "smarts": "[NH2][CH2][CH2][CH]([NH2])[C](=[O])[OH]"},
    {"name": "NMe-Ala", "smarts": "[CH3][NH][CH]([CH3])[C](=[O])[OH]"},
    {"name": "C9:0", "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][C](=[O])[OH]"},
    {
        "name": "NAc-COOH-Trp",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][c]1[c]([C](=[O])[OH])[nH][c]2[cH][cH][cH][cH][c]12)[C](=[O])[OH]",
    },
    {"name": "Ara/Lyx", "smarts": "[OH][CH]1[CH2][O][CH]([OH])[CH]([OH])[CH]1[OH]"},
    {
        "name": "NAc-Arg",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][CH2][CH2][N]=[C]([NH2])[NH2])[C](=[O])[OH]",
    },
    {
        "name": "NAc-C10:0-NH2(2)-oxo(8)",
        "smarts": "[CH3][CH2][C](=[O])[CH2][CH2][CH2][CH2][CH2][CH]([NH][C]([CH3])=[O])[C](=[O])[OH]",
    },
    {
        "name": "NMe-OMe-Ile",
        "smarts": "[CH3][CH2][CH]([CH3])[CH]([NH][CH3])[C](=[O])[O][CH3]",
    },
    {
        "name": "bOH-His",
        "smarts": "[NH2][C]([OH])([CH2][C](=[O])[OH])[c]1[cH][n][cH][nH]1",
    },
    {
        "name": "aC13:2(2.t4)",
        "smarts": "[CH3][CH2][CH]([CH3])[CH2][CH2][CH2][CH2][CH]=[CH][CH]=[CH][C](=[O])[OH]",
    },
    {
        "name": "NMe-Cl-OH-Trp",
        "smarts": "[CH3][NH][CH]([CH2][c]1[cH][nH][c]2[cH][c]([Cl])[c]([OH])[cH][c]12)[C](=[O])[OH]",
    },
    {"name": "bDpr", "smarts": "[NH2][CH]([NH2])[CH2][C](=[O])[OH]"},
    {
        "name": "C9:0-OH(3)",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH]([OH])[CH2][C](=[O])[OH]",
    },
    {"name": "bbOMe-Asp", "smarts": "[CH3][O][C](=[O])[CH]([NH2])[CH2][C](=[O])[OH]"},
    {
        "name": "NFo-Daz",
        "smarts": "[NH2][CH]([CH2][CH2][CH2][CH]([NH][CH]=[O])[C](=[O])[OH])[CH]([OH])[CH2][C](=[O])[OH]",
    },
    {
        "name": "F-ph-Gly",
        "smarts": "[NH2][CH]([C](=[O])[OH])[c]1[cH][cH][c]([F])[cH][cH]1",
    },
    {
        "name": "NAc-O-Met",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][CH2][S]([CH3])=[O])[C](=[O])[OH]",
    },
    {"name": "Oli", "smarts": "[CH3][CH]1[O][CH]([OH])[CH2][CH]([OH])[CH]1[OH]"},
    {
        "name": "ChrAct",
        "smarts": "[CH3][c]1[c]2[o][c]3[c]([CH3])[cH][cH][c]([C](=[O])[OH])[c]3[n][c]-2[c]([C](=[O])[OH])[c]([NH2])[c]1=[O]",
    },
    {
        "name": "NAc-Lys",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][CH2][CH2][CH2][NH2])[C](=[O])[OH]",
    },
    {
        "name": "NAc-Ahv",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][CH2][CH2][c]1[cH][cH][c]([OH])[cH][cH]1)[C](=[O])[OH]",
    },
    {
        "name": "Aca",
        "smarts": "[NH2][CH]([CH2][CH]1[CH2][CH2][C](=[O])[CH]2[O][CH]12)[C](=[O])[OH]",
    },
    {"name": "baFo-Gly", "smarts": "[NH2][C](=[O])[CH2][C](=[O])[OH]"},
    {
        "name": "NAc-4oxo-Hpr",
        "smarts": "[CH3][C](=[O])[N]1[CH2][CH2][C](=[O])[CH2][CH]1[C](=[O])[OH]",
    },
    {"name": "C4:0-OH(2.3.4)", "smarts": "[O]=[C]([OH])[CH]([OH])[CH]([OH])[CH2][OH]"},
    {
        "name": "bbOH-Gln",
        "smarts": "[NH2][C](=[O])[CH2][C]([NH2])([OH])[CH2][C](=[O])[OH]",
    },
    {"name": "NMe-Leu", "smarts": "[CH3][NH][CH]([CH2][CH]([CH3])[CH3])[C](=[O])[OH]"},
    {
        "name": "NMe-Aca",
        "smarts": "[CH3][NH][CH]([CH2][CH]1[CH2][CH2][C](=[O])[CH]2[O][CH]12)[C](=[O])[OH]",
    },
    {
        "name": "dh-Trp",
        "smarts": "[NH2][C](=[CH][c]1[cH][nH][c]2[cH][cH][cH][cH][c]12)[C](=[O])[OH]",
    },
    {
        "name": "NMe-Kyn",
        "smarts": "[CH3][NH][CH]([CH2][C](=[O])[c]1[cH][cH][cH][cH][c]1[NH2])[C](=[O])[OH]",
    },
    {
        "name": "Cl3-2OH-NMe-Leu",
        "smarts": "[CH3][NH][C]([OH])([CH2][CH]([CH3])[C]([Cl])([Cl])[Cl])[C](=[O])[OH]",
    },
    {
        "name": "Doe",
        "smarts": "[NH2][CH]([CH2][c]1[cH][cH][cH][cH][cH]1)[c]1[n][cH][cH][s]1",
    },
    {
        "name": "NAc-Har",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][CH2][CH2][CH2][N]=[C]([NH2])[NH2])[C](=[O])[OH]",
    },
    {"name": "bbMe-Asp", "smarts": "[CH3][C]([NH2])([CH2][C](=[O])[OH])[C](=[O])[OH]"},
    {
        "name": "bOH-NMe-Val",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[C]([CH3])([CH3])[OH]",
    },
    {
        "name": "NMe-Lys",
        "smarts": "[CH3][NH][CH]([CH2][CH2][CH2][CH2][NH2])[C](=[O])[OH]",
    },
    {
        "name": "C8:1(7)-OH(2.4.5)-NH2(3)-Ph(8)",
        "smarts": "[NH2][CH]([CH]([OH])[C](=[O])[OH])[CH]([OH])[CH]([OH])[CH2][CH]=[CH][c]1[cH][cH][cH][cH][cH]1",
    },
    {
        "name": "NAc-Cap",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[CH]1[CH2][CH2][N]=[C]([NH2])[NH]1",
    },
    {
        "name": "NOMe-Ac-Phe",
        "smarts": "[CH3][O][CH2][C](=[O])[NH][CH]([CH2][c]1[cH][cH][cH][cH][cH]1)[C](=[O])[OH]",
    },
    {
        "name": "dv-Tyr",
        "smarts": "[NH2][C]([CH]=[CH][C](=[O])[OH])=[CH][c]1[cH][cH][c]([OH])[cH][cH]1",
    },
    {
        "name": "ck-Arg",
        "smarts": "[NH]=[C]1[NH][C](=[O])[C]2([OH])[CH]([NH2])[CH2][CH2][CH2][N]12",
    },
    {
        "name": "NMe-Br-Trp",
        "smarts": "[CH3][NH][CH]([CH2][c]1[cH][nH][c]2[cH][cH][c]([Br])[cH][c]12)[C](=[O])[OH]",
    },
    {"name": "Gly", "smarts": "[NH2][CH2][C](=[O])[OH]"},
    {"name": "NMe-Orn", "smarts": "[CH3][NH][CH]([CH2][CH2][CH2][NH2])[C](=[O])[OH]"},
    {
        "name": "NMe-Dhpg",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[c]1[cH][c]([OH])[cH][c]([OH])[cH]1",
    },
    {
        "name": "NMe-Cl-Ile",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[CH]([CH3])[CH]([CH3])[Cl]",
    },
    {
        "name": "hk-Arg",
        "smarts": "[NH]=[C]([NH2])[NH][CH2][CH2][CH2][CH]([NH2])[CH]([OH])[C](=[O])[OH]",
    },
    {"name": "Me-Vaa", "smarts": "[CH3][CH2][CH]([CH3])[CH2][C](=[O])[OH]"},
    {
        "name": "NFo-pTrp",
        "smarts": "[O]=[CH][N]1[CH]([C](=[O])[OH])[CH2][C]2([OH])[c]3[cH][cH][cH][cH][c]3[NH][CH]12",
    },
    {
        "name": "ChrD",
        "smarts": "[NH2][CH]1[CH2][c]2[cH][c]([OH])[c]([OH])[cH][c]2[N]2[CH]([C](=[O])[OH])[CH2][CH2][NH][CH]12",
    },
    {"name": "Map", "smarts": "[CH3][CH2][CH]([NH2])[CH]([CH3])[C](=[O])[OH]"},
    {
        "name": "NFo-Ahv",
        "smarts": "[O]=[CH][NH][CH]([CH2][CH2][CH2][c]1[cH][cH][c]([OH])[cH][cH]1)[C](=[O])[OH]",
    },
    {"name": "Pda", "smarts": "[O]=[C]([OH])[CH2][CH2][CH2][C](=[O])[OH]"},
    {
        "name": "Hph",
        "smarts": "[NH2][CH]([CH2][CH2][c]1[cH][cH][cH][cH][cH]1)[C](=[O])[OH]",
    },
    {
        "name": "NAc-Ile",
        "smarts": "[CH3][CH2][CH]([CH3])[CH]([NH][C]([CH3])=[O])[C](=[O])[OH]",
    },
    {
        "name": "2OMe-Rha",
        "smarts": "[CH3][O][CH]1[CH]([OH])[O][CH]([CH3])[CH]([OH])[CH]1[OH]",
    },
    {"name": "bOH-Orn", "smarts": "[NH2][CH]([CH2][CH2][NH][OH])[CH2][C](=[O])[OH]"},
    {
        "name": "C8:0:1(7)-Me(2)",
        "smarts": "[CH]#[C][CH2][CH2][CH2][CH2][CH]([CH3])[C](=[O])[OH]",
    },
    {"name": "Nva", "smarts": "[CH3][CH2][CH2][CH]([NH2])[C](=[O])[OH]"},
    {
        "name": "NAc-Br-Tyr",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][c]1[cH][cH][c]([OH])[c]([Br])[cH]1)[C](=[O])[OH]",
    },
    {
        "name": "NFo-3OH-5Me-Pro",
        "smarts": "[CH3][CH]1[CH2][CH]([OH])[CH]([C](=[O])[OH])[N]1[CH]=[O]",
    },
    {
        "name": "aC13:0-OH(3)",
        "smarts": "[CH3][CH2][CH]([CH3])[CH2][CH2][CH2][CH2][CH2][CH2][CH]([OH])[CH2][C](=[O])[OH]",
    },
    {
        "name": "NAc-bbMe2-O-Met",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[C]([CH3])([CH3])[CH2][S]([CH3])=[O]",
    },
    {
        "name": "NMe-Br-OH-Trp",
        "smarts": "[CH3][NH][CH]([CH2][c]1[c]([Br])[nH][c]2[cH][cH][c]([OH])[cH][c]12)[C](=[O])[OH]",
    },
    {
        "name": "C6:0-OMe(3)",
        "smarts": "[CH3][CH2][CH2][CH]([CH2][C](=[O])[OH])[O][CH3]",
    },
    {
        "name": "iC12:0-OH(3)",
        "smarts": "[CH3][CH]([CH3])[CH2][CH2][CH2][CH2][CH2][CH2][CH]([OH])[CH2][C](=[O])[OH]",
    },
    {"name": "b4OH-Thr", "smarts": "[NH2][C]([OH])([CH2][OH])[CH2][C](=[O])[OH]"},
    {"name": "bIle", "smarts": "[CH3][CH2][C]([CH3])([NH2])[CH2][C](=[O])[OH]"},
    {
        "name": "NMe-bbMe2-O-Met",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[C]([CH3])([CH3])[CH2][S]([CH3])=[O]",
    },
    {
        "name": "DMOG",
        "smarts": "[CH3][CH]1[O][C]([c]2[cH][cH][cH][c]([OH])[c]2[OH])=[N][CH]1[C](=[O])[OH]",
    },
    {"name": "Ac-Aib", "smarts": "[CH3][C](=[O])[NH][C]([CH3])([CH3])[C](=[O])[OH]"},
    {
        "name": "NAc-3Me-4Me-Gln",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[CH]([CH3])[CH]([CH3])[C]([NH2])=[O]",
    },
    {
        "name": "bPO-Asn",
        "smarts": "[NH2][C](=[O])[C]([NH2])([CH2][C](=[O])[OH])[O][P](=[O])([OH])[OH]",
    },
    {
        "name": "NFo-OH-His",
        "smarts": "[O]=[CH][NH][CH]([C](=[O])[OH])[CH]([OH])[c]1[cH][n][cH][nH]1",
    },
    {
        "name": "NFo-Cl-Ile",
        "smarts": "[CH3][CH]([Cl])[CH]([CH3])[CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "NMe-diOH-Phe",
        "smarts": "[CH3][NH][CH]([CH2][c]1[cH][cH][c]([OH])[c]([OH])[cH]1)[C](=[O])[OH]",
    },
    {
        "name": "C8:3(t2.t4.t6)",
        "smarts": "[CH3][CH]=[CH][CH]=[CH][CH]=[CH][C](=[O])[OH]",
    },
    {
        "name": "NMe-OMe-Tyr",
        "smarts": "[CH3][NH][CH]([CH2][c]1[cH][cH][c]([O][CH3])[cH][cH]1)[C](=[O])[OH]",
    },
    {
        "name": "C8:0-OH(3)",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH]([OH])[CH2][C](=[O])[OH]",
    },
    {
        "name": "NFo-4OH-Thr",
        "smarts": "[O]=[CH][NH][CH]([C](=[O])[OH])[CH]([OH])[CH2][OH]",
    },
    {
        "name": "C9:1(Me4)-Me(2.4.6)-OH(8)-Oxo(5)",
        "smarts": "[CH2]=[C]([CH2][CH]([CH3])[C](=[O])[OH])[C](=[O])[CH]([CH3])[CH2][CH]([CH3])[OH]",
    },
    {
        "name": "NFo-diOH-Phe",
        "smarts": "[O]=[CH][NH][CH]([CH2][c]1[cH][cH][c]([OH])[c]([OH])[cH]1)[C](=[O])[OH]",
    },
    {
        "name": "NFo-PO-Asn",
        "smarts": "[NH2][C](=[O])[CH]([O][P](=[O])([OH])[OH])[CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "NdMe-Leu",
        "smarts": "[CH3][CH]([CH3])[CH2][CH]([C](=[O])[OH])[N]([CH3])[CH3]",
    },
    {
        "name": "NFo-Br-Tyr",
        "smarts": "[O]=[CH][NH][CH]([CH2][c]1[cH][cH][c]([OH])[c]([Br])[cH]1)[C](=[O])[OH]",
    },
    {
        "name": "NMe-C10:0-NH2(2)-oxo(8)",
        "smarts": "[CH3][CH2][C](=[O])[CH2][CH2][CH2][CH2][CH2][CH]([NH][CH3])[C](=[O])[OH]",
    },
    {
        "name": "C9:1(8)-Me(2)",
        "smarts": "[CH2]=[CH][CH2][CH2][CH2][CH2][CH2][CH]([CH3])[C](=[O])[OH]",
    },
    {
        "name": "NFo-N2Me-bOH-Asn",
        "smarts": "[CH3][NH][C](=[O])[CH]([OH])[CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "NFo-OH-Asn",
        "smarts": "[NH2][C](=[O])[CH]([OH])[CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "bbOH-Tyr",
        "smarts": "[NH2][C]([OH])([CH2][C](=[O])[OH])[c]1[cH][cH][c]([OH])[cH][cH]1",
    },
    {
        "name": "C10:0-OH(3)",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH]([OH])[CH2][C](=[O])[OH]",
    },
    {"name": "OH-Asn", "smarts": "[NH2][C](=[O])[CH]([OH])[CH]([NH2])[C](=[O])[OH]"},
    {
        "name": "iC11:0",
        "smarts": "[CH3][CH]([CH3])[CH2][CH2][CH2][CH2][CH2][CH2][CH2][C](=[O])[OH]",
    },
    {
        "name": "bN2Me-bOH-Asn",
        "smarts": "[CH3][NH][C](=[O])[C]([NH2])([OH])[CH2][C](=[O])[OH]",
    },
    {"name": "Ser", "smarts": "[NH2][CH]([CH2][OH])[C](=[O])[OH]"},
    {
        "name": "C10:0-NH2(2)-oxo(8)",
        "smarts": "[CH3][CH2][C](=[O])[CH2][CH2][CH2][CH2][CH2][CH]([NH2])[C](=[O])[OH]",
    },
    {
        "name": "C15:0-OH(3)",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH]([OH])[CH2][C](=[O])[OH]",
    },
    {
        "name": "NFo-3Me-Pro",
        "smarts": "[CH3][CH]1[CH2][CH2][N]([CH]=[O])[CH]1[C](=[O])[OH]",
    },
    {"name": "N2Me-Asn", "smarts": "[CH3][NH][C](=[O])[CH2][CH]([NH2])[C](=[O])[OH]"},
    {
        "name": "NAc-Hpg",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[c]1[cH][cH][c]([OH])[cH][cH]1",
    },
    {
        "name": "NAc-4Cl-Thr",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[CH]([OH])[CH2][Cl]",
    },
    {
        "name": "bEnd",
        "smarts": "[NH]=[C]1[NH][CH2][CH]([CH]([NH2])[CH2][C](=[O])[OH])[NH]1",
    },
    {
        "name": "Ac-Ival",
        "smarts": "[CH3][CH2][C]([CH3])([NH][C]([CH3])=[O])[C](=[O])[OH]",
    },
    {
        "name": "bOH-Trp",
        "smarts": "[NH2][CH]([CH2][C](=[O])[OH])[c]1[cH][nH][c]2[cH][cH][c]([OH])[cH][c]12",
    },
    {
        "name": "C13:0-NH2(3)",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH]([NH2])[CH2][C](=[O])[OH]",
    },
    {
        "name": "C8:0:1(7)-Me(2)-NH2(3)",
        "smarts": "[CH]#[C][CH2][CH2][CH2][CH]([NH2])[CH]([CH3])[C](=[O])[OH]",
    },
    {"name": "OH-cOrn", "smarts": "[NH2][CH]1[CH2][CH2][CH2][N]([OH])[C]1=[O]"},
    {
        "name": "bCl-CONH2-Trp",
        "smarts": "[NH2][C](=[O])[n]1[cH][c]([CH]([NH2])[CH2][C](=[O])[OH])[c]2[cH][cH][c]([Cl])[cH][c]21",
    },
    {
        "name": "C10:0-Me(2)-OH(3)",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH]([OH])[CH]([CH3])[C](=[O])[OH]",
    },
    {
        "name": "NMe-bMe-Gln",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[CH]([CH3])[CH2][C]([NH2])=[O]",
    },
    {"name": "bThr", "smarts": "[CH3][C]([NH2])([OH])[CH2][C](=[O])[OH]"},
    {"name": "C6:0-OH(3)", "smarts": "[CH3][CH2][CH2][CH]([OH])[CH2][C](=[O])[OH]"},
    {"name": "NFo-Asp", "smarts": "[O]=[CH][NH][CH]([CH2][C](=[O])[OH])[C](=[O])[OH]"},
    {
        "name": "NMe-bMe-Phe",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[CH]([CH3])[c]1[cH][cH][cH][cH][cH]1",
    },
    {
        "name": "C10:0:1(9)-Me(2.4)",
        "smarts": "[CH]#[C][CH2][CH2][CH2][CH2][CH]([CH3])[CH2][CH]([CH3])[C](=[O])[OH]",
    },
    {
        "name": "C12:1(5)",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH]=[CH][CH2][CH2][CH2][C](=[O])[OH]",
    },
    {
        "name": "NMe-Ahv",
        "smarts": "[CH3][NH][CH]([CH2][CH2][CH2][c]1[cH][cH][c]([OH])[cH][cH]1)[C](=[O])[OH]",
    },
    {"name": "Suc", "smarts": "[O]=[C]([OH])[CH2][CH2][C](=[O])[OH]"},
    {
        "name": "NFo-Cap",
        "smarts": "[NH2][C]1=[N][CH2][CH2][CH]([CH]([NH][CH]=[O])[C](=[O])[OH])[NH]1",
    },
    {
        "name": "iC12:2(2.t4)",
        "smarts": "[CH3][CH]([CH3])[CH2][CH2][CH2][CH2][CH]=[CH][CH]=[CH][C](=[O])[OH]",
    },
    {"name": "NFo-Val", "smarts": "[CH3][CH]([CH3])[CH]([NH][CH]=[O])[C](=[O])[OH]"},
    {
        "name": "bMe-Br-Phe",
        "smarts": "[CH3][CH]([c]1[cH][cH][c]([Br])[cH][cH]1)[CH]([NH2])[C](=[O])[OH]",
    },
    {"name": "Lac", "smarts": "[CH3][CH]([OH])[C](=[O])[OH]"},
    {
        "name": "NAc-diOH-Arg",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[CH]([OH])[CH]([OH])[CH2][N]=[C]([NH2])[NH2]",
    },
    {
        "name": "NAc-C10:0-NH2(2)-Ep(9)-oxo(8)",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][CH2][CH2][CH2][CH2][C](=[O])[CH]1[CH2][O]1)[C](=[O])[OH]",
    },
    {
        "name": "NMe-Cl2-Hpg",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[c]1[cH][c]([Cl])[c]([OH])[c]([Cl])[cH]1",
    },
    {
        "name": "I-NMe-Tyr",
        "smarts": "[CH3][NH][CH]([CH2][c]1[cH][cH][c]([OH])[c]([I])[cH]1)[C](=[O])[OH]",
    },
    {
        "name": "C18:1(9)-OH(3)",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH]=[CH][CH2][CH2][CH2][CH2][CH2][CH]([OH])[CH2][C](=[O])[OH]",
    },
    {
        "name": "Hip",
        "smarts": "[CH3][CH]([C](=[O])[OH])[C](=[O])[CH]([OH])[CH]([CH3])[CH3]",
    },
    {
        "name": "NFo-COOH-Trp",
        "smarts": "[O]=[CH][NH][CH]([CH2][c]1[c]([C](=[O])[OH])[nH][c]2[cH][cH][cH][cH][c]12)[C](=[O])[OH]",
    },
    {"name": "pGlu", "smarts": "[O]=[C]1[CH2][CH2][CH]([C](=[O])[OH])[NH]1"},
    {
        "name": "PO-Asn",
        "smarts": "[NH2][C](=[O])[CH]([O][P](=[O])([OH])[OH])[CH]([NH2])[C](=[O])[OH]",
    },
    {
        "name": "NAc-5OH-Cap",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[CH]1[CH2][CH]([OH])[NH][C](=[NH])[NH]1",
    },
    {
        "name": "NMe-Ahad",
        "smarts": "[CH3][NH][CH]([CH2][CH]([OH])[CH2][C](=[O])[OH])[C](=[O])[OH]",
    },
    {
        "name": "C14:0",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][C](=[O])[OH]",
    },
    {
        "name": "NMe-3NO2-Tyr",
        "smarts": "[CH3][NH][CH]([CH2][c]1[cH][cH][c]([OH])[c]([N+](=[O])[O-])[cH]1)[C](=[O])[OH]",
    },
    {
        "name": "NMe-MeO-Glu",
        "smarts": "[CH3][NH][CH]([CH2][CH2][C](=[O])[O][CH3])[C](=[O])[OH]",
    },
    {"name": "bOMe-Thr", "smarts": "[CH3][O][C]([CH3])([NH2])[CH2][C](=[O])[OH]"},
    {
        "name": "bBr-OH-Trp",
        "smarts": "[NH2][CH]([CH2][C](=[O])[OH])[c]1[c]([Br])[nH][c]2[cH][cH][c]([OH])[cH][c]12",
    },
    {
        "name": "Cl2-Hpg",
        "smarts": "[NH2][CH]([C](=[O])[OH])[c]1[cH][c]([Cl])[c]([OH])[c]([Cl])[cH]1",
    },
    {"name": "Asp", "smarts": "[NH2][CH]([CH2][C](=[O])[OH])[C](=[O])[OH]"},
    {"name": "L-Dbu/Dbu", "smarts": "[CH3][CH]([NH2])[CH]([NH2])[C](=[O])[OH]"},
    {
        "name": "NFo-Trp",
        "smarts": "[O]=[CH][NH][CH]([CH2][c]1[cH][nH][c]2[cH][cH][cH][cH][c]12)[C](=[O])[OH]",
    },
    {
        "name": "ChrA",
        "smarts": "[O]=[C]([OH])[CH]1[CH2][CH2][N]2[C](=[O])[NH][C]3=[CH][c]4[cH][c]([OH])[c]([OH])[cH][c]4[N]1[CH]32",
    },
    {
        "name": "NAc-bOH-Br-Phe",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[CH]([OH])[c]1[cH][cH][c]([Br])[cH][cH]1",
    },
    {
        "name": "iC10:0",
        "smarts": "[CH3][CH]([CH3])[CH2][CH2][CH2][CH2][CH2][CH2][C](=[O])[OH]",
    },
    {
        "name": "NAc-Ac-OH-Orn",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][CH2][CH2][N]([OH])[C]([CH3])=[O])[C](=[O])[OH]",
    },
    {
        "name": "aC10:0",
        "smarts": "[CH3][CH2][CH]([CH3])[CH2][CH2][CH2][CH2][CH2][C](=[O])[OH]",
    },
    {
        "name": "NMe-3Me-Phe",
        "smarts": "[CH3][NH][CH]([CH2][c]1[cH][cH][cH][c]([CH3])[cH]1)[C](=[O])[OH]",
    },
    {
        "name": "NFo-Ph-Gly",
        "smarts": "[O]=[CH][NH][CH]([C](=[O])[OH])[c]1[cH][cH][cH][cH][cH]1",
    },
    {
        "name": "bCit",
        "smarts": "[NH2][C](=[O])[NH][CH2][CH2][CH]([NH2])[CH2][C](=[O])[OH]",
    },
    {"name": "2Dh-Mabu", "smarts": "[CH3][CH]=[C]([NH][CH3])[C](=[O])[OH]"},
    {"name": "Hpoe", "smarts": "[O]=[C]([OH])[C](=[O])[c]1[cH][cH][c]([OH])[cH][cH]1"},
    {"name": "4OH-Thr", "smarts": "[NH2][CH]([C](=[O])[OH])[CH]([OH])[CH2][OH]"},
    {
        "name": "NMe-bMe-Leu/diMe-aIle",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[CH]([CH3])[CH]([CH3])[CH3]",
    },
    {
        "name": "DMAdda",
        "smarts": "[CH3][C]([CH]=[CH][CH]([NH2])[CH]([CH3])[C](=[O])[OH])=[CH][CH]([CH3])[CH]([OH])[CH2][c]1[cH][cH][cH][cH][cH]1",
    },
    {
        "name": "NAc-Dab",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][CH2][NH2])[C](=[O])[OH]",
    },
    {"name": "diMe-Cys", "smarts": "[CH3][NH][CH]([CH2][S][CH3])[C](=[O])[OH]"},
    {"name": "C8:2(2.t4)", "smarts": "[CH3][CH2][CH2][CH]=[CH][CH]=[CH][C](=[O])[OH]"},
    {
        "name": "NFo-Br-OH-Trp",
        "smarts": "[O]=[CH][NH][CH]([CH2][c]1[c]([Br])[nH][c]2[cH][cH][c]([OH])[cH][c]12)[C](=[O])[OH]",
    },
    {
        "name": "NFo-Amv",
        "smarts": "[CH3][O][c]1[cH][cH][c]([CH2][CH2][CH2][CH]([NH][CH]=[O])[C](=[O])[OH])[cH][cH]1",
    },
    {
        "name": "NAc-3NO2-Tyr",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][c]1[cH][cH][c]([OH])[c]([N+](=[O])[O-])[cH]1)[C](=[O])[OH]",
    },
    {
        "name": "iC9:0-OH(3)",
        "smarts": "[CH3][CH]([CH3])[CH2][CH2][CH2][CH]([OH])[CH2][C](=[O])[OH]",
    },
    {
        "name": "NAc-4oxo-Pro",
        "smarts": "[CH3][C](=[O])[N]1[CH2][C](=[O])[CH2][CH]1[C](=[O])[OH]",
    },
    {
        "name": "NMe-Et-Nva",
        "smarts": "[CH3][CH2][CH]([CH2][CH3])[CH]([NH][CH3])[C](=[O])[OH]",
    },
    {"name": "bOMe-Asp", "smarts": "[CH3][O][C](=[O])[CH2][CH]([NH2])[C](=[O])[OH]"},
    {"name": "bLeu", "smarts": "[CH3][CH]([CH3])[CH]([NH2])[CH2][C](=[O])[OH]"},
    {"name": "Cl-Ile", "smarts": "[CH3][CH]([Cl])[CH]([CH3])[CH]([NH2])[C](=[O])[OH]"},
    {
        "name": "bNMe-Lan",
        "smarts": "[NH2][CH]([CH2][C](=[O])[OH])[S][CH2][CH]([NH][OH])[C](=[O])[OH]",
    },
    {
        "name": "NFo-bMe-Br-Phe",
        "smarts": "[CH3][CH]([c]1[cH][cH][c]([Br])[cH][cH]1)[CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "NAc-Hty",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][CH2][c]1[cH][cH][c]([OH])[cH][cH]1)[C](=[O])[OH]",
    },
    {
        "name": "NMe-Phe",
        "smarts": "[CH3][NH][CH]([CH2][c]1[cH][cH][cH][cH][cH]1)[C](=[O])[OH]",
    },
    {
        "name": "NMe-diOH-Arg",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[CH]([OH])[CH]([OH])[CH2][N]=[C]([NH2])[NH2]",
    },
    {
        "name": "NAc-diOH-Phe",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][c]1[cH][cH][c]([OH])[c]([OH])[cH]1)[C](=[O])[OH]",
    },
    {
        "name": "bBr-Phe",
        "smarts": "[NH2][CH]([CH2][C](=[O])[OH])[c]1[cH][cH][c]([Br])[cH][cH]1",
    },
    {
        "name": "Cl-Trp",
        "smarts": "[NH2][CH]([CH2][C]1=[N][c]2[cH][c]([Cl])[cH][cH][c]2[CH2]1)[C](=[O])[OH]",
    },
    {"name": "Pha", "smarts": "[O]=[C]([OH])[CH2][c]1[cH][cH][cH][cH][cH]1"},
    {"name": "NFo-Hse", "smarts": "[O]=[CH][NH][CH]([CH2][CH2][OH])[C](=[O])[OH]"},
    {
        "name": "NMe-C10:0-OH(8)-NH2(2)",
        "smarts": "[CH3][CH2][CH]([OH])[CH2][CH2][CH2][CH2][CH2][CH]([NH][CH3])[C](=[O])[OH]",
    },
    {
        "name": "bAmv",
        "smarts": "[CH3][O][c]1[cH][cH][c]([CH2][CH2][CH]([NH2])[CH2][C](=[O])[OH])[cH][cH]1",
    },
    {
        "name": "bOMe-Trp",
        "smarts": "[CH3][O][n]1[cH][c]([CH]([NH2])[CH2][C](=[O])[OH])[c]2[cH][cH][cH][cH][c]21",
    },
    {
        "name": "C9:1(4)-Me(2.4.6)-OH(8)",
        "smarts": "[CH3][C](=[CH][CH]([CH3])[CH2][CH]([CH3])[OH])[CH2][CH]([CH3])[C](=[O])[OH]",
    },
    {
        "name": "iC10:2(2.t4)",
        "smarts": "[CH3][CH]([CH3])[CH2][CH2][CH]=[CH][CH]=[CH][C](=[O])[OH]",
    },
    {
        "name": "C13:0-OH(3)",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH]([OH])[CH2][C](=[O])[OH]",
    },
    {"name": "bMe-Asp", "smarts": "[CH3][CH]([C](=[O])[OH])[CH]([NH2])[C](=[O])[OH]"},
    {
        "name": "iC5:0-OH(2)-CA(4)",
        "smarts": "[CH3][CH]([C](=[O])[OH])[CH]([OH])[C](=[O])[OH]",
    },
    {
        "name": "NAc-bMe-Ile",
        "smarts": "[CH3][CH2][C]([CH3])([CH3])[CH]([NH][C]([CH3])=[O])[C](=[O])[OH]",
    },
    {"name": "Met", "smarts": "[CH3][S][CH2][CH2][CH]([NH2])[C](=[O])[OH]"},
    {
        "name": "iC16:0-NH2(3)",
        "smarts": "[CH3][CH]([CH3])[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH]([NH2])[CH2][C](=[O])[OH]",
    },
    {
        "name": "bdiOH-Arg",
        "smarts": "[NH2][C]([NH2])=[N][CH2][CH]([OH])[C]([NH2])([OH])[CH2][C](=[O])[OH]",
    },
    {
        "name": "NAc-CysA",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][S](=[O])(=[O])[OH])[C](=[O])[OH]",
    },
    {"name": "NMe-Gly-Thz", "smarts": "[CH3][NH][CH2][c]1[n][cH][cH][s]1"},
    {
        "name": "iC14:1(3)",
        "smarts": "[CH3][CH]([CH3])[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH]=[CH][CH2][C](=[O])[OH]",
    },
    {
        "name": "v-OH-Tyr",
        "smarts": "[NH2][CH]([CH]=[CH][C](=[O])[OH])[CH2][c]1[cH][cH][c]([OH])[c]([OH])[cH]1",
    },
    {
        "name": "NMe-4Me-Pro",
        "smarts": "[CH3][CH]1[CH2][CH]([C](=[O])[OH])[N]([CH3])[CH2]1",
    },
    {"name": "dh-Ala", "smarts": "[CH2]=[C]([NH2])[C](=[O])[OH]"},
    {"name": "4Me-Hva", "smarts": "[CH3][CH]([CH3])[CH2][CH]([OH])[C](=[O])[OH]"},
    {
        "name": "NFo-3NO2-Tyr",
        "smarts": "[O]=[CH][NH][CH]([CH2][c]1[cH][cH][c]([OH])[c]([N+](=[O])[O-])[cH]1)[C](=[O])[OH]",
    },
    {
        "name": "aC11:0",
        "smarts": "[CH3][CH2][CH]([CH3])[CH2][CH2][CH2][CH2][CH2][CH2][C](=[O])[OH]",
    },
    {
        "name": "NAc-Met",
        "smarts": "[CH3][S][CH2][CH2][CH]([NH][C]([CH3])=[O])[C](=[O])[OH]",
    },
    {"name": "NMe-Dbu", "smarts": "[CH3][NH][CH]([C](=[O])[OH])[CH]([CH3])[NH2]"},
    {"name": "Pro-Thz", "smarts": "[cH]1[cH][s][c]([CH]2[CH2][CH2][CH2][NH]2)[n]1"},
    {
        "name": "NAc-pGlu",
        "smarts": "[CH3][C](=[O])[N]1[C](=[O])[CH2][CH2][CH]1[C](=[O])[OH]",
    },
    {"name": "5Me-Pro", "smarts": "[CH3][CH]1[CH2][CH2][CH]([C](=[O])[OH])[NH]1"},
    {
        "name": "N2Me-bOH-Asn",
        "smarts": "[CH3][NH][C](=[O])[CH]([OH])[CH]([NH2])[C](=[O])[OH]",
    },
    {
        "name": "NMe-His",
        "smarts": "[CH3][NH][CH]([CH2][c]1[cH][n][cH][nH]1)[C](=[O])[OH]",
    },
    {
        "name": "NAc-End",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][CH]1[CH2][NH][C](=[NH])[NH]1)[C](=[O])[OH]",
    },
    {
        "name": "NAc-OH-Orn",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][CH2][CH2][NH][OH])[C](=[O])[OH]",
    },
    {"name": "dhAbu", "smarts": "[CH3][CH]=[C]([NH2])[C](=[O])[OH]"},
    {
        "name": "COOH-Trp",
        "smarts": "[NH2][CH]([CH2][c]1[c]([C](=[O])[OH])[nH][c]2[cH][cH][cH][cH][c]12)[C](=[O])[OH]",
    },
    {
        "name": "iC12:0",
        "smarts": "[CH3][CH]([CH3])[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][C](=[O])[OH]",
    },
    {
        "name": "NFo-bOH-Gln",
        "smarts": "[NH2][C](=[O])[CH2][CH]([OH])[CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {"name": "NAc-Cys", "smarts": "[CH3][C](=[O])[NH][CH]([CH2][SH])[C](=[O])[OH]"},
    {
        "name": "NFo-3OH-Leu",
        "smarts": "[CH3][CH]([CH3])[CH]([OH])[CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "NMe-Ac-Ser",
        "smarts": "[CH3][NH][CH]([CH2][O][C]([CH3])=[O])[C](=[O])[OH]",
    },
    {
        "name": "aC13:0",
        "smarts": "[CH3][CH2][CH]([CH3])[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][C](=[O])[OH]",
    },
    {
        "name": "bArg",
        "smarts": "[NH2][C]([NH2])=[N][CH2][CH2][CH]([NH2])[CH2][C](=[O])[OH]",
    },
    {
        "name": "C14:1(7)",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH]=[CH][CH2][CH2][CH2][CH2][CH2][C](=[O])[OH]",
    },
    {
        "name": "NAc-3OH-Pro",
        "smarts": "[CH3][C](=[O])[N]1[CH2][CH2][CH]([OH])[CH]1[C](=[O])[OH]",
    },
    {
        "name": "C15:0-NH2(3)",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH]([NH2])[CH2][C](=[O])[OH]",
    },
    {"name": "NMe-Dab", "smarts": "[CH3][NH][CH]([CH2][CH2][NH2])[C](=[O])[OH]"},
    {
        "name": "dDap",
        "smarts": "[CH3][CH]([C](=[O])[OH])[CH]([OH])[CH]1[CH2][CH2][CH2][NH]1",
    },
    {
        "name": "NFo-5Me-Pro",
        "smarts": "[CH3][CH]1[CH2][CH2][CH]([C](=[O])[OH])[N]1[CH]=[O]",
    },
    {
        "name": "NFo-Me-AOA",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "bAc-OH-Orn",
        "smarts": "[CH3][C](=[O])[N]([OH])[CH2][CH2][CH]([NH2])[CH2][C](=[O])[OH]",
    },
    {
        "name": "bMe-Phe",
        "smarts": "[CH3][CH]([c]1[cH][cH][cH][cH][cH]1)[CH]([NH2])[C](=[O])[OH]",
    },
    {
        "name": "b3NO2-Tyr",
        "smarts": "[NH2][CH]([CH2][C](=[O])[OH])[c]1[cH][cH][c]([OH])[c]([N+](=[O])[O-])[cH]1",
    },
    {
        "name": "Cl2-NMe-Leu",
        "smarts": "[CH3][NH][CH]([CH2][CH]([CH3])[CH]([Cl])[Cl])[C](=[O])[OH]",
    },
    {
        "name": "NMe-OH-Trp",
        "smarts": "[CH3][NH][CH]([CH2][c]1[cH][nH][c]2[cH][cH][c]([OH])[cH][c]12)[C](=[O])[OH]",
    },
    {
        "name": "bPhe",
        "smarts": "[NH2][CH]([CH2][C](=[O])[OH])[c]1[cH][cH][cH][cH][cH]1",
    },
    {"name": "NFo-Hpr", "smarts": "[O]=[CH][N]1[CH2][CH2][CH2][CH2][CH]1[C](=[O])[OH]"},
    {
        "name": "NAc-pTrp",
        "smarts": "[CH3][C](=[O])[N]1[CH]([C](=[O])[OH])[CH2][C]2([OH])[c]3[cH][cH][cH][cH][c]3[NH][CH]12",
    },
    {
        "name": "NAc-ChrI",
        "smarts": "[CH3][C](=[O])[N]1[CH]([C](=[O])[OH])[CH2][CH2][N]2[c]3[cH][c]([OH])[c]([OH])[cH][c]3[CH]=[C]([NH2])[CH]21",
    },
    {
        "name": "C10:0-Me(2.2.4)-OH(3.7)",
        "smarts": "[CH3][CH2][CH2][CH]([OH])[CH2][CH2][CH]([CH3])[CH]([OH])[C]([CH3])([CH3])[C](=[O])[OH]",
    },
    {"name": "NFo-Pro", "smarts": "[O]=[CH][N]1[CH2][CH2][CH2][CH]1[C](=[O])[OH]"},
    {"name": "gSer", "smarts": "[NH]=[C]([NH2])[NH][C]([NH2])([OH])[C](=[O])[OH]"},
    {
        "name": "b3OH-Leu",
        "smarts": "[CH3][CH]([CH3])[C]([NH2])([OH])[CH2][C](=[O])[OH]",
    },
    {
        "name": "iC8:0-Me(2.4)-OH(3)",
        "smarts": "[CH3][CH]([CH3])[CH2][CH]([CH3])[CH]([OH])[CH]([CH3])[C](=[O])[OH]",
    },
    {
        "name": "aC11:2(4.6)-Me(2.6)-OH(2.3)",
        "smarts": "[CH3][CH2][CH]([CH3])[CH]=[C]([CH3])[CH]=[CH][CH]([OH])[C]([CH3])([OH])[C](=[O])[OH]",
    },
    {
        "name": "NFo-N2Me-Asn",
        "smarts": "[CH3][NH][C](=[O])[CH2][CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "NFo-C10:0-OH(8)-NH2(2)",
        "smarts": "[CH3][CH2][CH]([OH])[CH2][CH2][CH2][CH2][CH2][CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "bMeO-Glu",
        "smarts": "[CH3][O][C](=[O])[CH2][CH]([NH2])[CH2][C](=[O])[OH]",
    },
    {
        "name": "NAc-Cl-Tyr",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][c]1[cH][cH][c]([OH])[c]([Cl])[cH]1)[C](=[O])[OH]",
    },
    {
        "name": "NMe-PT",
        "smarts": "[CH3][NH][CH]([CH2][CH2][P]([CH3])(=[O])[OH])[C](=[O])[OH]",
    },
    {
        "name": "NAc-Fo-OH-Orn",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][CH2][CH2][N]([OH])[CH]=[O])[C](=[O])[OH]",
    },
    {"name": "Gln", "smarts": "[NH2][C](=[O])[CH2][CH2][CH]([NH2])[C](=[O])[OH]"},
    {
        "name": "NAc-Orn",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][CH2][CH2][NH2])[C](=[O])[OH]",
    },
    {"name": "NAc-Abu", "smarts": "[CH3][CH2][CH]([NH][C]([CH3])=[O])[C](=[O])[OH]"},
    {
        "name": "NMe-Daz",
        "smarts": "[CH3][NH][CH]([CH2][CH2][CH2][CH]([NH2])[CH]([OH])[CH2][C](=[O])[OH])[C](=[O])[OH]",
    },
    {
        "name": "NMe-F-ph-Gly",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[c]1[cH][cH][c]([F])[cH][cH]1",
    },
    {
        "name": "NMe-bOH-Br-Phe",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[CH]([OH])[c]1[cH][cH][c]([Br])[cH][cH]1",
    },
    {
        "name": "bBr-Trp",
        "smarts": "[NH2][CH]([CH2][C](=[O])[OH])[c]1[cH][nH][c]2[cH][cH][c]([Br])[cH][c]12",
    },
    {"name": "Hpr", "smarts": "[O]=[C]([OH])[CH]1[CH2][CH2][CH2][CH2][NH]1"},
    {"name": "Me-Suc", "smarts": "[CH3][CH]1[CH2][C](=[O])[NH][C]1=[O]"},
    {
        "name": "NMe-4oxo-Pro",
        "smarts": "[CH3][N]1[CH2][C](=[O])[CH2][CH]1[C](=[O])[OH]",
    },
    {"name": "bLys", "smarts": "[NH2][CH2][CH2][CH2][CH]([NH2])[CH2][C](=[O])[OH]"},
    {
        "name": "NMe-Cl-Hpg",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[c]1[cH][cH][c]([OH])[c]([Cl])[cH]1",
    },
    {
        "name": "Br-NMe-Tyr",
        "smarts": "[CH3][NH][CH]([CH2][c]1[cH][cH][c]([OH])[c]([Br])[cH]1)[C](=[O])[OH]",
    },
    {
        "name": "C16:0",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][C](=[O])[OH]",
    },
    {"name": "iC7:0", "smarts": "[CH3][CH]([CH3])[CH2][CH2][CH2][C](=[O])[OH]"},
    {"name": "bDbu", "smarts": "[CH3][C]([NH2])([NH2])[CH2][C](=[O])[OH]"},
    {
        "name": "C12:1(11)-Me(6)-OH(2.4.5)-NH2(3)-mPhe(11)",
        "smarts": "[CH3][O][c]1[cH][cH][c]([CH]=[CH][CH2][CH2][CH2][CH2][CH]([CH3])[CH]([OH])[CH]([OH])[CH]([NH2])[CH]([OH])[C](=[O])[OH])[cH][cH]1",
    },
    {
        "name": "NMe-Br-Trp",
        "smarts": "[CH3][NH][CH]([CH2][c]1[c]([Br])[nH][c]2[cH][cH][cH][cH][c]12)[C](=[O])[OH]",
    },
    {"name": "NFo-Cys", "smarts": "[O]=[CH][NH][CH]([CH2][SH])[C](=[O])[OH]"},
    {
        "name": "C11:0-OH(3)",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH]([OH])[CH2][C](=[O])[OH]",
    },
    {
        "name": "NMe-2Me-3Me-pGlu",
        "smarts": "[CH3][CH]1[C](=[O])[N]([CH3])[CH]([C](=[O])[OH])[CH]1[CH3]",
    },
    {
        "name": "NFo-4OH-Pro",
        "smarts": "[O]=[CH][N]1[CH2][CH]([OH])[CH2][CH]1[C](=[O])[OH]",
    },
    {
        "name": "diOH-Phe",
        "smarts": "[NH2][CH]([CH2][c]1[cH][cH][c]([OH])[c]([OH])[cH]1)[C](=[O])[OH]",
    },
    {
        "name": "bbMe-Br-Phe",
        "smarts": "[CH3][C]([NH2])([CH2][C](=[O])[OH])[c]1[cH][cH][c]([Br])[cH][cH]1",
    },
    {
        "name": "bbMe-Phe",
        "smarts": "[CH3][C]([NH2])([CH2][C](=[O])[OH])[c]1[cH][cH][cH][cH][cH]1",
    },
    {"name": "diOH-Bz", "smarts": "[O]=[C]([OH])[c]1[cH][cH][cH][c]([OH])[c]1[OH]"},
    {"name": "Ileol", "smarts": "[CH3][CH2][CH]([CH3])[CH]([NH2])[CH2][OH]"},
    {
        "name": "NFo-3OH-Pro",
        "smarts": "[O]=[CH][N]1[CH2][CH2][CH]([OH])[CH]1[C](=[O])[OH]",
    },
    {"name": "Ile/aIle", "smarts": "[CH3][CH2][CH]([CH3])[CH]([NH2])[C](=[O])[OH]"},
    {"name": "bAsn", "smarts": "[NH2][C](=[O])[CH]([NH2])[CH2][C](=[O])[OH]"},
    {"name": "N-OH-Hta", "smarts": "[OH][NH][CH2][CH2][CH]1[CH2][NH][CH2][NH]1"},
    {
        "name": "Van/Ere",
        "smarts": "[CH3][CH]1[O][CH]([OH])[CH2][C]([CH3])([NH2])[CH]1[OH]",
    },
    {"name": "NAc-Ala", "smarts": "[CH3][C](=[O])[NH][CH]([CH3])[C](=[O])[OH]"},
    {
        "name": "bbOMe-Tyr",
        "smarts": "[CH3][O][C]([NH2])([CH2][C](=[O])[OH])[c]1[cH][cH][c]([OH])[cH][cH]1",
    },
    {
        "name": "C8:0-Me(2.2)-OH(3)",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH]([OH])[C]([CH3])([CH3])[C](=[O])[OH]",
    },
    {
        "name": "NMe-OH-Asn",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[CH]([OH])[C]([NH2])=[O]",
    },
    {
        "name": "NFo-OH-Orn",
        "smarts": "[O]=[CH][NH][CH]([CH2][CH2][CH2][NH][OH])[C](=[O])[OH]",
    },
    {"name": "C5:0-NH2(3)", "smarts": "[CH3][CH2][CH]([NH2])[CH2][C](=[O])[OH]"},
    {
        "name": "OH-Trp",
        "smarts": "[NH2][CH]([CH2][c]1[cH][nH][c]2[cH][cH][c]([OH])[cH][c]12)[C](=[O])[OH]",
    },
    {
        "name": "NMe-1Me-Trp",
        "smarts": "[CH3][NH][CH]([CH2][c]1[cH][n]([CH3])[c]2[cH][cH][cH][cH][c]12)[C](=[O])[OH]",
    },
    {
        "name": "bbMe-Gln",
        "smarts": "[CH3][C]([NH2])([CH2][C]([NH2])=[O])[CH2][C](=[O])[OH]",
    },
    {
        "name": "C15:2(t4.t6)-OH(2.3)",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH]=[CH][CH]=[CH][CH]([OH])[CH]([OH])[C](=[O])[OH]",
    },
    {"name": "C7:0", "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][C](=[O])[OH]"},
    {"name": "OH-Asp", "smarts": "[NH2][CH]([C](=[O])[OH])[CH]([OH])[C](=[O])[OH]"},
    {"name": "Et-Nva", "smarts": "[CH3][CH2][CH]([CH2][CH3])[CH]([NH2])[C](=[O])[OH]"},
    {"name": "NMe-Nva", "smarts": "[CH3][CH2][CH2][CH]([NH][CH3])[C](=[O])[OH]"},
    {
        "name": "NAc-Hph",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][CH2][c]1[cH][cH][cH][cH][cH]1)[C](=[O])[OH]",
    },
    {
        "name": "NMe-Tyr",
        "smarts": "[CH3][NH][CH]([CH2][c]1[cH][cH][c]([OH])[cH][cH]1)[C](=[O])[OH]",
    },
    {
        "name": "NAc-MeO-Glu",
        "smarts": "[CH3][O][C](=[O])[CH2][CH2][CH]([NH][C]([CH3])=[O])[C](=[O])[OH]",
    },
    {"name": "3OH-Leu", "smarts": "[CH3][CH]([CH3])[CH]([OH])[CH]([NH2])[C](=[O])[OH]"},
    {"name": "bAad", "smarts": "[NH2][CH]([CH2][CH2][C](=[O])[OH])[CH2][C](=[O])[OH]"},
    {
        "name": "NMe-NMe-Lan",
        "smarts": "[CH3][NH][CH]([CH2][S][CH2][CH]([NH][OH])[C](=[O])[OH])[C](=[O])[OH]",
    },
    {"name": "Aad", "smarts": "[NH2][CH]([CH2][CH2][CH2][C](=[O])[OH])[C](=[O])[OH]"},
    {"name": "CO", "smarts": "[CH2]=[O]"},
    {"name": "NFo-Abu", "smarts": "[CH3][CH2][CH]([NH][CH]=[O])[C](=[O])[OH]"},
    {
        "name": "NMe-Cl-CONH2-Trp",
        "smarts": "[CH3][NH][CH]([CH2][c]1[cH][n]([C]([NH2])=[O])[c]2[cH][c]([Cl])[cH][cH][c]12)[C](=[O])[OH]",
    },
    {
        "name": "NMe-O-Met",
        "smarts": "[CH3][NH][CH]([CH2][CH2][S]([CH3])=[O])[C](=[O])[OH]",
    },
    {
        "name": "Dhpg",
        "smarts": "[NH2][CH]([C](=[O])[OH])[c]1[cH][c]([OH])[cH][c]([OH])[cH]1",
    },
    {
        "name": "bCOOH-Trp",
        "smarts": "[NH2][CH]([CH2][C](=[O])[OH])[c]1[c]([C](=[O])[OH])[nH][c]2[cH][cH][cH][cH][c]12",
    },
    {
        "name": "bCl-Tyr",
        "smarts": "[NH2][CH]([CH2][C](=[O])[OH])[c]1[cH][cH][c]([OH])[c]([Cl])[cH]1",
    },
    {
        "name": "NMe-ChrI",
        "smarts": "[CH3][N]1[CH]([C](=[O])[OH])[CH2][CH2][N]2[c]3[cH][c]([OH])[c]([OH])[cH][c]3[CH]=[C]([NH2])[CH]21",
    },
    {"name": "OMe-Thr", "smarts": "[CH3][O][CH]([CH3])[CH]([NH2])[C](=[O])[OH]"},
    {
        "name": "bOH-Cl-Tyr",
        "smarts": "[NH2][CH]([C](=[O])[OH])[CH]([OH])[c]1[cH][cH][c]([OH])[c]([Cl])[cH]1",
    },
    {"name": "bAc-Ser", "smarts": "[CH3][C](=[O])[O][CH]([NH2])[CH2][C](=[O])[OH]"},
    {
        "name": "ph-Gly/Ph-Gly",
        "smarts": "[NH2][CH]([C](=[O])[OH])[c]1[cH][cH][cH][cH][cH]1",
    },
    {"name": "Ac-Pro", "smarts": "[CH3][C](=[O])[N]1[CH2][CH2][CH2][CH]1[C](=[O])[OH]"},
    {
        "name": "NFo-Hpg",
        "smarts": "[O]=[CH][NH][CH]([C](=[O])[OH])[c]1[cH][cH][c]([OH])[cH][cH]1",
    },
    {
        "name": "bdiOH-Phe",
        "smarts": "[NH2][CH]([CH2][C](=[O])[OH])[c]1[cH][cH][c]([OH])[c]([OH])[cH]1",
    },
    {
        "name": "NAc-Hse",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][CH2][OH])[C](=[O])[OH]",
    },
    {
        "name": "aC15:0",
        "smarts": "[CH3][CH2][CH]([CH3])[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][C](=[O])[OH]",
    },
    {
        "name": "C7:0-Me(2)-OH(3)",
        "smarts": "[CH3][CH2][CH2][CH2][CH]([OH])[CH]([CH3])[C](=[O])[OH]",
    },
    {
        "name": "NAc-Kyn",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][C](=[O])[c]1[cH][cH][cH][cH][c]1[NH2])[C](=[O])[OH]",
    },
    {
        "name": "NFo-5OH-Cap",
        "smarts": "[NH]=[C]1[NH][CH]([OH])[CH2][CH]([CH]([NH][CH]=[O])[C](=[O])[OH])[NH]1",
    },
    {
        "name": "C16:1(7)",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH]=[CH][CH2][CH2][CH2][CH2][CH2][C](=[O])[OH]",
    },
    {
        "name": "GSpd",
        "smarts": "[NH2][CH2][CH2][CH2][NH][CH2][CH2][CH2][CH2][N]=[C]([NH2])[NH2]",
    },
    {
        "name": "Me2-Bmt",
        "smarts": "[CH3][CH]=[CH][CH2][C]([CH3])([CH3])[CH]([OH])[CH]([NH][CH3])[C](=[O])[OH]",
    },
    {
        "name": "NFo-Choi",
        "smarts": "[O]=[CH][N]1[CH]([C](=[O])[OH])[CH2][CH]2[CH2][CH2][CH]([OH])[CH2][CH]21",
    },
    {
        "name": "bBr-Tyr",
        "smarts": "[NH2][CH]([CH2][C](=[O])[OH])[c]1[cH][cH][c]([OH])[c]([Br])[cH]1",
    },
    {"name": "CysA", "smarts": "[NH2][CH]([CH2][S](=[O])(=[O])[OH])[C](=[O])[OH]"},
    {
        "name": "Br-Trp",
        "smarts": "[NH2][CH]([CH2][c]1[cH][nH][c]2[cH][cH][c]([Br])[cH][c]12)[C](=[O])[OH]",
    },
    {
        "name": "NMe-OMe-Trp",
        "smarts": "[CH3][NH][CH]([CH2][c]1[cH][n]([O][CH3])[c]2[cH][cH][cH][cH][c]12)[C](=[O])[OH]",
    },
    {"name": "NMe-Asn", "smarts": "[CH3][NH][CH]([CH2][C]([NH2])=[O])[C](=[O])[OH]"},
    {
        "name": "Cl3-5OH-NMe-Leu",
        "smarts": "[CH3][NH][CH]([CH2][CH]([CH2][OH])[C]([Cl])([Cl])[Cl])[C](=[O])[OH]",
    },
    {
        "name": "NFo-Ac-OH-Orn",
        "smarts": "[CH3][C](=[O])[N]([OH])[CH2][CH2][CH2][CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {"name": "C4:0-OH(2.3)", "smarts": "[CH3][CH]([OH])[CH]([OH])[C](=[O])[OH]"},
    {
        "name": "NFo-OH-Trp",
        "smarts": "[O]=[CH][NH][CH]([CH2][c]1[cH][nH][c]2[cH][cH][c]([OH])[cH][c]12)[C](=[O])[OH]",
    },
    {
        "name": "iC13:0-OH(3)",
        "smarts": "[CH3][CH]([CH3])[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH]([OH])[CH2][C](=[O])[OH]",
    },
    {
        "name": "bOH-Gln",
        "smarts": "[NH2][C](=[O])[CH2][CH]([OH])[CH]([NH2])[C](=[O])[OH]",
    },
    {
        "name": "Fo-OH-Orn",
        "smarts": "[NH2][CH]([CH2][CH2][CH2][N]([OH])[CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "aC15:1(3)",
        "smarts": "[CH3][CH2][CH]([CH3])[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH]=[CH][CH2][C](=[O])[OH]",
    },
    {
        "name": "C6:0-Me(2)-NH2(3)",
        "smarts": "[CH3][CH2][CH2][CH]([NH2])[CH]([CH3])[C](=[O])[OH]",
    },
    {
        "name": "iC17:0-NH2(3)",
        "smarts": "[CH3][CH]([CH3])[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH]([NH2])[CH2][C](=[O])[OH]",
    },
    {
        "name": "COOH-Qui",
        "smarts": "[O]=[C]([OH])[c]1[cH][n][c]2[cH][cH][cH][cH][c]2[n]1",
    },
    {"name": "NMe-Hse", "smarts": "[CH3][NH][CH]([CH2][CH2][OH])[C](=[O])[OH]"},
    {"name": "NFo-aFo-Gly", "smarts": "[O]=[CH][NH][CH]([CH]=[O])[C](=[O])[OH]"},
    {
        "name": "iC13:0",
        "smarts": "[CH3][CH]([CH3])[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][C](=[O])[OH]",
    },
    {
        "name": "NAc-Ph-Gly",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[c]1[cH][cH][cH][cH][cH]1",
    },
    {
        "name": "NAc-Cit",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][CH2][CH2][NH][C]([NH2])=[O])[C](=[O])[OH]",
    },
    {
        "name": "Hysp",
        "smarts": "[CH3][CH2][CH]([CH3])[CH]([OH])[C](=[O])[CH]([CH3])[C](=[O])[OH]",
    },
    {
        "name": "NFo-3Me-4Me-Gln",
        "smarts": "[CH3][CH]([C]([NH2])=[O])[CH]([CH3])[CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {"name": "Act", "smarts": "[CH3][O][CH]1[CH]([NH2])[CH2][CH]([OH])[O][CH]1[CH3]"},
    {"name": "NMe-bAla", "smarts": "[CH3][NH][CH2][CH2][C](=[O])[OH]"},
    {"name": "bHil", "smarts": "[CH3][CH2][CH]([CH3])[CH]([NH2])[CH2][C](=[O])[OH]"},
    {
        "name": "NFo-C10:0-NH2(2)-Ep(9)-oxo(8)",
        "smarts": "[O]=[CH][NH][CH]([CH2][CH2][CH2][CH2][CH2][C](=[O])[CH]1[CH2][O]1)[C](=[O])[OH]",
    },
    {"name": "Tig", "smarts": "[CH3][CH]=[C]([CH3])[C](=[O])[OH]"},
    {"name": "His", "smarts": "[NH2][CH]([CH2][c]1[cH][n][cH][nH]1)[C](=[O])[OH]"},
    {
        "name": "Bmt",
        "smarts": "[CH3][CH]=[CH][CH2][CH]([CH3])[CH]([OH])[CH]([NH2])[C](=[O])[OH]",
    },
    {"name": "Ival", "smarts": "[CH3][CH2][C]([CH3])([NH2])[C](=[O])[OH]"},
    {"name": "4OH-Pro", "smarts": "[O]=[C]([OH])[CH]1[CH2][CH]([OH])[CH2][NH]1"},
    {
        "name": "Arg",
        "smarts": "[NH2][C]([NH2])=[N][CH2][CH2][CH2][CH]([NH2])[C](=[O])[OH]",
    },
    {"name": "Put", "smarts": "[NH2][CH2][CH2][CH2][CH2][NH2]"},
    {
        "name": "NMe-Aad",
        "smarts": "[CH3][NH][CH]([CH2][CH2][CH2][C](=[O])[OH])[C](=[O])[OH]",
    },
    {"name": "3Me-Pro", "smarts": "[CH3][CH]1[CH2][CH2][NH][CH]1[C](=[O])[OH]"},
    {
        "name": "NAc-Et-Nva",
        "smarts": "[CH3][CH2][CH]([CH2][CH3])[CH]([NH][C]([CH3])=[O])[C](=[O])[OH]",
    },
    {"name": "NFo-Dbu", "smarts": "[CH3][CH]([NH2])[CH]([NH][CH]=[O])[C](=[O])[OH]"},
    {
        "name": "PALOA",
        "smarts": "[CH3][CH]([NH2])[c]1[n][c]([CH]=[CH][C](=[O])[OH])[cH][o]1",
    },
    {
        "name": "NAc-OH-Trp",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][c]1[cH][nH][c]2[cH][cH][c]([OH])[cH][c]12)[C](=[O])[OH]",
    },
    {
        "name": "DHPT",
        "smarts": "[O]=[C]([OH])[CH]1[CH2][S][C]([c]2[cH][cH][cH][c]([OH])[c]2[OH])=[N]1",
    },
    {"name": "Iser", "smarts": "[NH2][CH2][CH]([OH])[C](=[O])[OH]"},
    {
        "name": "C6:0-Me(2.2)-OH(3)",
        "smarts": "[CH3][CH2][CH2][CH]([OH])[C]([CH3])([CH3])[C](=[O])[OH]",
    },
    {
        "name": "bbMe-NMe-Trp",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[C]([CH3])([CH3])[c]1[cH][nH][c]2[cH][cH][cH][cH][c]12",
    },
    {"name": "CMA", "smarts": "[CH3][CH2][CH]1[CH2][C]1([NH2])[C](=[O])[OH]"},
    {
        "name": "NFo-Cl2-Hpg",
        "smarts": "[O]=[CH][NH][CH]([C](=[O])[OH])[c]1[cH][c]([Cl])[c]([OH])[c]([Cl])[cH]1",
    },
    {
        "name": "NMe-N2Me-bOH-Asn",
        "smarts": "[CH3][NH][C](=[O])[CH]([OH])[CH]([NH][CH3])[C](=[O])[OH]",
    },
    {
        "name": "bAhv",
        "smarts": "[NH2][CH]([CH2][CH2][c]1[cH][cH][c]([OH])[cH][cH]1)[CH2][C](=[O])[OH]",
    },
    {
        "name": "NMe-MeA-Phe",
        "smarts": "[CH3][NH][c]1[cH][cH][c]([CH2][CH]([NH][CH3])[C](=[O])[OH])[cH][cH]1",
    },
    {
        "name": "b3Me-Glu",
        "smarts": "[CH3][C]([NH2])([CH2][C](=[O])[OH])[CH2][C](=[O])[OH]",
    },
    {"name": "pNH2-Bz", "smarts": "[NH2][c]1[cH][cH][c]([C](=[O])[OH])[cH][cH]1"},
    {
        "name": "dDil",
        "smarts": "[CH3][CH2][CH]([CH3])[CH]([NH2])[CH]([CH2][C](=[O])[OH])[O][CH3]",
    },
    {
        "name": "NFo-Fo-OH-Orn",
        "smarts": "[O]=[CH][NH][CH]([CH2][CH2][CH2][N]([OH])[CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "NMe-Cap",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[CH]1[CH2][CH2][N]=[C]([NH2])[NH]1",
    },
    {
        "name": "Cit",
        "smarts": "[NH2][C](=[O])[NH][CH2][CH2][CH2][CH]([NH2])[C](=[O])[OH]",
    },
    {
        "name": "NAc-Aad",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][CH2][CH2][C](=[O])[OH])[C](=[O])[OH]",
    },
    {
        "name": "NAc-OMe-Trp",
        "smarts": "[CH3][O][n]1[cH][c]([CH2][CH]([NH][C]([CH3])=[O])[C](=[O])[OH])[c]2[cH][cH][cH][cH][c]21",
    },
    {
        "name": "NAc-bMe-Asp",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[CH]([CH3])[C](=[O])[OH]",
    },
    {
        "name": "Cl-CONH2-Trp",
        "smarts": "[NH2][C](=[O])[n]1[cH][c]([CH2][CH]([NH2])[C](=[O])[OH])[c]2[cH][cH][c]([Cl])[cH][c]21",
    },
    {
        "name": "NMe-3OH-Pro",
        "smarts": "[CH3][N]1[CH2][CH2][CH]([OH])[CH]1[C](=[O])[OH]",
    },
    {
        "name": "OAc-Leuol",
        "smarts": "[CH3][C](=[O])[O][CH2][CH]([NH2])[CH2][CH]([CH3])[CH3]",
    },
    {"name": "bGln", "smarts": "[NH2][C](=[O])[CH2][CH]([NH2])[CH2][C](=[O])[OH]"},
    {"name": "Orn", "smarts": "[NH2][CH2][CH2][CH2][CH]([NH2])[C](=[O])[OH]"},
    {
        "name": "Dpy",
        "smarts": "[CH3][O][C]1=[CH][C](=[O])[NH][CH]1[CH2][c]1[cH][cH][cH][cH][cH]1",
    },
    {
        "name": "3Me-Hty",
        "smarts": "[CH3][CH]([CH2][c]1[cH][cH][c]([OH])[cH][cH]1)[CH]([NH2])[C](=[O])[OH]",
    },
    {
        "name": "NFo-diOH-Arg",
        "smarts": "[NH2][C]([NH2])=[N][CH2][CH]([OH])[CH]([OH])[CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {"name": "bVal", "smarts": "[CH3][C]([CH3])([NH2])[CH2][C](=[O])[OH]"},
    {
        "name": "NFo-Ahad",
        "smarts": "[O]=[CH][NH][CH]([CH2][CH]([OH])[CH2][C](=[O])[OH])[C](=[O])[OH]",
    },
    {
        "name": "Cap",
        "smarts": "[NH2][C]1=[N][CH2][CH2][CH]([CH]([NH2])[C](=[O])[OH])[NH]1",
    },
    {"name": "NFo-Ser", "smarts": "[O]=[CH][NH][CH]([CH2][OH])[C](=[O])[OH]"},
    {"name": "NMe-Val", "smarts": "[CH3][NH][CH]([C](=[O])[OH])[CH]([CH3])[CH3]"},
    {
        "name": "NMe-Choi",
        "smarts": "[CH3][N]1[CH]([C](=[O])[OH])[CH2][CH]2[CH2][CH2][CH]([OH])[CH2][CH]21",
    },
    {
        "name": "iC9:0",
        "smarts": "[CH3][CH]([CH3])[CH2][CH2][CH2][CH2][CH2][C](=[O])[OH]",
    },
    {
        "name": "Ac-OH-Orn",
        "smarts": "[CH3][C](=[O])[N]([OH])[CH2][CH2][CH2][CH]([NH2])[C](=[O])[OH]",
    },
    {
        "name": "C8:2(5.7)-Me(6)-OH(4)-NH2(3)-brPh(8)",
        "smarts": "[CH3][C]([CH]=[CH][c]1[cH][cH][c]([Br])[cH][cH]1)=[CH][CH]([OH])[CH]([NH2])[CH2][C](=[O])[OH]",
    },
    {
        "name": "4oxo-Van",
        "smarts": "[CH3][CH]1[O][CH]([OH])[CH2][C]([CH3])([NH2])[C]1=[O]",
    },
    {
        "name": "Ac-Phe",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][c]1[cH][cH][cH][cH][cH]1)[C](=[O])[OH]",
    },
    {
        "name": "NMe-bMe-Ile",
        "smarts": "[CH3][CH2][C]([CH3])([CH3])[CH]([NH][CH3])[C](=[O])[OH]",
    },
    {
        "name": "bMe-Ile",
        "smarts": "[CH3][CH2][C]([CH3])([CH3])[CH]([NH2])[C](=[O])[OH]",
    },
    {
        "name": "C10:0-Me(2.2.4)-OH(3)-OMe(7)",
        "smarts": "[CH3][CH2][CH2][CH]([CH2][CH2][CH]([CH3])[CH]([OH])[C]([CH3])([CH3])[C](=[O])[OH])[O][CH3]",
    },
    {
        "name": "bCap",
        "smarts": "[NH2][C]1=[N][CH2][CH2][C]([NH2])([CH2][C](=[O])[OH])[NH]1",
    },
    {
        "name": "NFo-Hty",
        "smarts": "[O]=[CH][NH][CH]([CH2][CH2][c]1[cH][cH][c]([OH])[cH][cH]1)[C](=[O])[OH]",
    },
    {"name": "iC8:0", "smarts": "[CH3][CH]([CH3])[CH2][CH2][CH2][CH2][C](=[O])[OH]"},
    {"name": "MCP", "smarts": "[CH3][n]1[cH][c]([Cl])[cH][c]1[C](=[O])[OH]"},
    {
        "name": "C14:2(t4.t6)-OH(2.3)",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH]=[CH][CH]=[CH][CH]([OH])[CH]([OH])[C](=[O])[OH]",
    },
    {
        "name": "NAc-His",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][c]1[cH][n][cH][nH]1)[C](=[O])[OH]",
    },
    {"name": "C6:2(t2.t4)", "smarts": "[CH3][CH]=[CH][CH]=[CH][C](=[O])[OH]"},
    {
        "name": "NMe-Cit",
        "smarts": "[CH3][NH][CH]([CH2][CH2][CH2][NH][C]([NH2])=[O])[C](=[O])[OH]",
    },
    {
        "name": "NAc-Choi",
        "smarts": "[CH3][C](=[O])[N]1[CH]([C](=[O])[OH])[CH2][CH]2[CH2][CH2][CH]([OH])[CH2][CH]21",
    },
    {
        "name": "bKyn",
        "smarts": "[NH2][c]1[cH][cH][cH][cH][c]1[C](=[O])[CH]([NH2])[CH2][C](=[O])[OH]",
    },
    {
        "name": "NFo-Ile",
        "smarts": "[CH3][CH2][CH]([CH3])[CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "NFo-Tyr",
        "smarts": "[O]=[CH][NH][CH]([CH2][c]1[cH][cH][c]([OH])[cH][cH]1)[C](=[O])[OH]",
    },
    {
        "name": "NAc-Cl-Ile",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[CH]([CH3])[CH]([CH3])[Cl]",
    },
    {
        "name": "Ahp",
        "smarts": "[CH3][CH]([OH])[CH]([C](=[O])[OH])[N]1[C](=[O])[CH]([NH2])[CH2][CH2][CH]1[OH]",
    },
    {
        "name": "NMe-3Me-4Me-Gln",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[CH]([CH3])[CH]([CH3])[C]([NH2])=[O]",
    },
    {
        "name": "ChrI",
        "smarts": "[NH2][C]1=[CH][c]2[cH][c]([OH])[c]([OH])[cH][c]2[N]2[CH2][CH2][CH]([C](=[O])[OH])[NH][CH]12",
    },
    {"name": "bHse", "smarts": "[NH2][CH]([CH2][OH])[CH2][C](=[O])[OH]"},
    {
        "name": "aC17:0-NH2(3)",
        "smarts": "[CH3][CH2][CH]([CH3])[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH]([NH2])[CH2][C](=[O])[OH]",
    },
    {"name": "Mab", "smarts": "[CH3][CH]([NH2])[CH]([CH3])[C](=[O])[OH]"},
    {"name": "Hmp", "smarts": "[CH3][CH2][CH]([CH3])[CH]([OH])[C](=[O])[OH]"},
    {"name": "C4:1(3)-OH(2)", "smarts": "[CH2]=[CH][CH]([OH])[C](=[O])[OH]"},
    {"name": "Ivalol", "smarts": "[CH3][CH2][C]([CH3])([NH2])[CH2][OH]"},
    {
        "name": "NFo-Hil",
        "smarts": "[CH3][CH2][CH]([CH3])[CH2][CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "PTTA",
        "smarts": "[NH2][CH]([CH2][c]1[cH][cH][c]([OH])[cH][cH]1)[c]1[n][c]([CH]=[CH][C](=[O])[OH])[cH][s]1",
    },
    {
        "name": "bC10:0-OH(9)-NH2(2)",
        "smarts": "[CH3][CH]([OH])[CH2][CH2][CH2][CH2][CH2][CH]([NH2])[CH2][C](=[O])[OH]",
    },
    {
        "name": "OH-dHpg",
        "smarts": "[O]=[C]([OH])[C](=[N][OH])[c]1[cH][cH][c]([OH])[cH][cH]1",
    },
    {
        "name": "v-Arg",
        "smarts": "[NH]=[C]([NH2])[NH][CH2][CH2][CH2][CH]([NH2])[CH]=[CH][C](=[O])[OH]",
    },
    {
        "name": "b1Me-Trp",
        "smarts": "[CH3][n]1[cH][c]([CH]([NH2])[CH2][C](=[O])[OH])[c]2[cH][cH][cH][cH][c]21",
    },
    {
        "name": "NMe-bOMe-Asp",
        "smarts": "[CH3][NH][CH]([CH2][C](=[O])[O][CH3])[C](=[O])[OH]",
    },
    {
        "name": "iC14:0-NH2(3)",
        "smarts": "[CH3][CH]([CH3])[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH]([NH2])[CH2][C](=[O])[OH]",
    },
    {"name": "Ria/Aco", "smarts": "[CH3][CH]1[O][CH]([OH])[CH2][CH]([NH2])[CH]1[OH]"},
    {"name": "Pheol", "smarts": "[NH2][CH]([CH2][OH])[CH2][c]1[cH][cH][cH][cH][cH]1"},
    {
        "name": "NAc-Cl-Hpg",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[c]1[cH][cH][c]([OH])[c]([Cl])[cH]1",
    },
    {
        "name": "NAc-Me-AOA",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH]([NH][C]([CH3])=[O])[C](=[O])[OH]",
    },
    {"name": "Ala-Thz", "smarts": "[CH3][CH]([NH2])[c]1[n][c]([C](=[O])[OH])[cH][s]1"},
    {
        "name": "NFo-End",
        "smarts": "[NH]=[C]1[NH][CH2][CH]([CH2][CH]([NH][CH]=[O])[C](=[O])[OH])[NH]1",
    },
    {
        "name": "Amv",
        "smarts": "[CH3][O][c]1[cH][cH][c]([CH2][CH2][CH2][CH]([NH2])[C](=[O])[OH])[cH][cH]1",
    },
    {
        "name": "NMe-bOH-Tyr",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[CH]([OH])[c]1[cH][cH][c]([OH])[cH][cH]1",
    },
    {
        "name": "NAc-C10:0-OH(9)-NH2(2)-oxo(8)",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][CH2][CH2][CH2][CH2][C](=[O])[CH]([CH3])[OH])[C](=[O])[OH]",
    },
    {
        "name": "NMe-N2Me-Asn",
        "smarts": "[CH3][NH][C](=[O])[CH2][CH]([NH][CH3])[C](=[O])[OH]",
    },
    {
        "name": "NAc-Cl2-Pro",
        "smarts": "[CH3][C](=[O])[N]1[CH2][CH]([Cl])[CH]([Cl])[CH]1[C](=[O])[OH]",
    },
    {"name": "t-Leu", "smarts": "[CH3][C]([CH3])([CH3])[CH]([NH2])[C](=[O])[OH]"},
    {
        "name": "Har",
        "smarts": "[N]#[C][N]([CH2][CH2][CH2][CH2][CH]([OH])[C]([NH2])=[O])[N]=[O]",
    },
    {"name": "pOH-Bz", "smarts": "[O]=[C]([OH])[c]1[cH][cH][c]([OH])[cH][cH]1"},
    {
        "name": "bApv",
        "smarts": "[NH2][CH]([CH2][CH2][c]1[cH][cH][cH][cH][cH]1)[CH2][C](=[O])[OH]",
    },
    {
        "name": "NFo-Cl-Tyr",
        "smarts": "[O]=[CH][NH][CH]([CH2][c]1[cH][cH][c]([OH])[c]([Cl])[cH]1)[C](=[O])[OH]",
    },
    {"name": "C4:0-OH(2)-Ep(3)", "smarts": "[O]=[C]([OH])[CH]([OH])[CH]1[CH2][O]1"},
    {
        "name": "NMe-Cl-Trp",
        "smarts": "[CH3][NH][CH]([CH2][c]1[cH][nH][c]2[cH][c]([Cl])[cH][cH][c]12)[C](=[O])[OH]",
    },
    {
        "name": "NAc-bOH-Gln",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[CH]([OH])[CH2][C]([NH2])=[O]",
    },
    {
        "name": "NMe-4OH-Thr",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[CH]([OH])[CH2][OH]",
    },
    {"name": "3OMe-Ala", "smarts": "[CH3][O][C](=[O])[CH]([CH3])[NH2]"},
    {
        "name": "NAc-Bmt",
        "smarts": "[CH3][CH]=[CH][CH2][CH]([CH3])[CH]([OH])[CH]([NH][C]([CH3])=[O])[C](=[O])[OH]",
    },
    {
        "name": "NMe-O2-Met",
        "smarts": "[CH3][NH][CH]([CH2][CH2][S]([CH3])(=[O])=[O])[C](=[O])[OH]",
    },
    {
        "name": "NAc-Leu",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][CH]([CH3])[CH3])[C](=[O])[OH]",
    },
    {
        "name": "3Me-4Me-Gln",
        "smarts": "[CH3][CH]([C]([NH2])=[O])[CH]([CH3])[CH]([NH2])[C](=[O])[OH]",
    },
    {
        "name": "NMe-3OH-Leu",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[CH]([OH])[CH]([CH3])[CH3]",
    },
    {
        "name": "NMe-bOH-Cl-Tyr",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[CH]([OH])[c]1[cH][cH][c]([OH])[c]([Cl])[cH]1",
    },
    {"name": "C4:0-Me(2)", "smarts": "[CH3][CH2][CH]([CH3])[C](=[O])[OH]"},
    {
        "name": "Cl-NMe-Tyr",
        "smarts": "[CH3][NH][CH]([CH2][c]1[cH][cH][c]([OH])[c]([Cl])[cH]1)[C](=[O])[OH]",
    },
    {
        "name": "NFo-Arg",
        "smarts": "[NH2][C]([NH2])=[N][CH2][CH2][CH2][CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "NMe-5OH-Cap",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[CH]1[CH2][CH]([OH])[NH][C](=[NH])[NH]1",
    },
    {
        "name": "NAc-3OH-Leu",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[CH]([OH])[CH]([CH3])[CH3]",
    },
    {"name": "Valol", "smarts": "[CH3][CH]([CH3])[CH]([NH2])[CH2][OH]"},
    {"name": "NAc-aFo-Gly", "smarts": "[CH3][C](=[O])[NH][CH]([CH]=[O])[C](=[O])[OH]"},
    {
        "name": "NMe-Amv",
        "smarts": "[CH3][NH][CH]([CH2][CH2][CH2][c]1[cH][cH][c]([O][CH3])[cH][cH]1)[C](=[O])[OH]",
    },
    {"name": "NMe-Dha", "smarts": "[CH2]=[C]([NH][CH3])[C](=[O])[OH]"},
    {
        "name": "C10:0-OH(8)-NH2(2)",
        "smarts": "[CH3][CH2][CH]([OH])[CH2][CH2][CH2][CH2][CH2][CH]([NH2])[C](=[O])[OH]",
    },
    {"name": "Leu", "smarts": "[CH3][CH]([CH3])[CH2][CH]([NH2])[C](=[O])[OH]"},
    {"name": "NFo-Met", "smarts": "[CH3][S][CH2][CH2][CH]([NH][CH]=[O])[C](=[O])[OH]"},
    {"name": "bO-Met", "smarts": "[CH3][S](=[O])[CH2][CH]([NH2])[CH2][C](=[O])[OH]"},
    {
        "name": "NAc-Amv",
        "smarts": "[CH3][O][c]1[cH][cH][c]([CH2][CH2][CH2][CH]([NH][C]([CH3])=[O])[C](=[O])[OH])[cH][cH]1",
    },
    {
        "name": "NOMe-Ac-Val",
        "smarts": "[CH3][O][CH2][C](=[O])[NH][CH]([C](=[O])[OH])[CH]([CH3])[CH3]",
    },
    {
        "name": "End",
        "smarts": "[NH]=[C]1[NH][CH2][CH]([CH2][CH]([NH2])[C](=[O])[OH])[NH]1",
    },
    {
        "name": "Br-Phe",
        "smarts": "[NH2][CH]([CH2][c]1[cH][cH][c]([Br])[cH][cH]1)[C](=[O])[OH]",
    },
    {"name": "Hse", "smarts": "[NH2][CH]([CH2][CH2][OH])[C](=[O])[OH]"},
    {
        "name": "NMe-Cl-Trp",
        "smarts": "[CH3][NH][CH]([CH2][C]1=[N][c]2[cH][c]([Cl])[cH][cH][c]2[CH2]1)[C](=[O])[OH]",
    },
    {"name": "bCysA", "smarts": "[NH2][CH]([CH2][C](=[O])[OH])[S](=[O])(=[O])[OH]"},
    {
        "name": "C14:0-OH(3.4)",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH]([OH])[CH]([OH])[CH2][C](=[O])[OH]",
    },
    {
        "name": "Trp",
        "smarts": "[NH2][CH]([CH2][c]1[cH][nH][c]2[cH][cH][cH][cH][c]12)[C](=[O])[OH]",
    },
    {"name": "OH-Pyr", "smarts": "[O]=[C]1[CH2][CH]([OH])[CH2][NH]1"},
    {
        "name": "Sta",
        "smarts": "[CH3][CH]([CH3])[CH2][CH]([NH2])[CH]([OH])[CH2][C](=[O])[OH]",
    },
    {
        "name": "3NO2-Tyr",
        "smarts": "[NH2][CH]([CH2][c]1[cH][cH][c]([OH])[c]([N+](=[O])[O-])[cH]1)[C](=[O])[OH]",
    },
    {
        "name": "NFo-Cl-Hpg",
        "smarts": "[O]=[CH][NH][CH]([C](=[O])[OH])[c]1[cH][cH][c]([OH])[c]([Cl])[cH]1",
    },
    {"name": "C8:0", "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][C](=[O])[OH]"},
    {
        "name": "NMe-bOMe-Tyr",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[CH]([O][CH3])[c]1[cH][cH][c]([OH])[cH][cH]1",
    },
    {
        "name": "NMe-Glu",
        "smarts": "[CH3][NH][CH]([CH2][CH2][C](=[O])[OH])[C](=[O])[OH]",
    },
    {
        "name": "NFo-4oxo-Pro",
        "smarts": "[O]=[CH][N]1[CH2][C](=[O])[CH2][CH]1[C](=[O])[OH]",
    },
    {
        "name": "NAc-5Me-Pro",
        "smarts": "[CH3][C](=[O])[N]1[CH]([CH3])[CH2][CH2][CH]1[C](=[O])[OH]",
    },
    {"name": "Pya", "smarts": "[CH3][C](=[O])[C](=[O])[OH]"},
    {
        "name": "NFo-O2-Met",
        "smarts": "[CH3][S](=[O])(=[O])[CH2][CH2][CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {"name": "bSer", "smarts": "[NH2][CH]([OH])[CH2][C](=[O])[OH]"},
    {"name": "bCys", "smarts": "[NH2][CH]([SH])[CH2][C](=[O])[OH]"},
    {
        "name": "C8:0:1(7)-Me(2.2)-OH(3)",
        "smarts": "[CH]#[C][CH2][CH2][CH2][CH]([OH])[C]([CH3])([CH3])[C](=[O])[OH]",
    },
    {
        "name": "PAOA",
        "smarts": "[CH3][CH2][CH]([NH2])[c]1[n][c]([CH]=[CH][C](=[O])[OH])[cH][o]1",
    },
    {"name": "NMe-aFo-Gly", "smarts": "[CH3][NH][CH]([CH]=[O])[C](=[O])[OH]"},
    {
        "name": "b3Me-4Me-Gln",
        "smarts": "[CH3][CH]([C]([NH2])=[O])[C]([CH3])([NH2])[CH2][C](=[O])[OH]",
    },
    {
        "name": "Cl3-NMe-dhLeu",
        "smarts": "[CH3][NH][C](=[CH][CH]([CH3])[C]([Cl])([Cl])[Cl])[C](=[O])[OH]",
    },
    {
        "name": "NFo-Leu",
        "smarts": "[CH3][CH]([CH3])[CH2][CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "NAc-C10:0-OH(9)-NH2(2)",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][CH2][CH2][CH2][CH2][CH2][CH]([CH3])[OH])[C](=[O])[OH]",
    },
    {"name": "Pro", "smarts": "[O]=[C]([OH])[CH]1[CH2][CH2][CH2][NH]1"},
    {
        "name": "NMe-4oxo-Hpr",
        "smarts": "[CH3][N]1[CH2][CH2][C](=[O])[CH2][CH]1[C](=[O])[OH]",
    },
    {
        "name": "bO2-Met",
        "smarts": "[CH3][S](=[O])(=[O])[CH2][CH]([NH2])[CH2][C](=[O])[OH]",
    },
    {
        "name": "PMST",
        "smarts": "[CH3][O][CH2][CH]([NH2])[c]1[n][c]([CH]=[CH][C](=[O])[OH])[cH][s]1",
    },
    {"name": "Leuol", "smarts": "[CH3][CH]([CH3])[CH2][CH]([NH2])[CH2][OH]"},
    {"name": "4Cl-Thr", "smarts": "[NH2][CH]([C](=[O])[OH])[CH]([OH])[CH2][Cl]"},
    {
        "name": "NMe-Ph-Gly",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[c]1[cH][cH][cH][cH][cH]1",
    },
    {
        "name": "C8:1(7)-Me(2.2)-OH(3)",
        "smarts": "[CH2]=[CH][CH2][CH2][CH2][CH]([OH])[C]([CH3])([CH3])[C](=[O])[OH]",
    },
    {"name": "NFo-Thr", "smarts": "[CH3][CH]([OH])[CH]([NH][CH]=[O])[C](=[O])[OH]"},
    {
        "name": "NMe-OH-Asp",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[CH]([OH])[C](=[O])[OH]",
    },
    {"name": "dhCys", "smarts": "[NH2][C](=[CH][SH])[C](=[O])[OH]"},
    {
        "name": "NAc-NMe-Lan",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][S][CH2][CH]([NH][OH])[C](=[O])[OH])[C](=[O])[OH]",
    },
    {
        "name": "NMe-Cl2-Pro",
        "smarts": "[CH3][N]1[CH2][CH]([Cl])[CH]([Cl])[CH]1[C](=[O])[OH]",
    },
    {"name": "Pyr", "smarts": "[O]=[C]1[CH2][CH2][CH2][NH]1"},
    {
        "name": "bTrp",
        "smarts": "[NH2][CH]([CH2][C](=[O])[OH])[c]1[cH][nH][c]2[cH][cH][cH][cH][c]12",
    },
    {"name": "4Me-Pro", "smarts": "[CH3][CH]1[CH2][NH][CH]([C](=[O])[OH])[CH2]1"},
    {
        "name": "NAc-N2Me-Asn",
        "smarts": "[CH3][NH][C](=[O])[CH2][CH]([NH][C]([CH3])=[O])[C](=[O])[OH]",
    },
    {
        "name": "NMe-3Me-Glu",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[CH]([CH3])[CH2][C](=[O])[OH]",
    },
    {
        "name": "Adda",
        "smarts": "[CH3][O][CH]([CH2][c]1[cH][cH][cH][cH][cH]1)[CH]([CH3])[CH]=[C]([CH3])[CH]=[CH][CH]([NH2])[CH]([CH3])[C](=[O])[OH]",
    },
    {
        "name": "NMe-Apv",
        "smarts": "[CH3][NH][CH]([CH2][CH2][CH2][c]1[cH][cH][cH][cH][cH]1)[C](=[O])[OH]",
    },
    {
        "name": "NMe-4oxo-5Me-Pro",
        "smarts": "[CH3][CH]1[C](=[O])[CH2][CH]([C](=[O])[OH])[N]1[CH3]",
    },
    {
        "name": "Ahv",
        "smarts": "[NH2][CH]([CH2][CH2][CH2][c]1[cH][cH][c]([OH])[cH][cH]1)[C](=[O])[OH]",
    },
    {
        "name": "NMe-dPhe",
        "smarts": "[CH3][NH][C](=[CH][c]1[cH][cH][cH][cH][cH]1)[C](=[O])[OH]",
    },
    {
        "name": "NMe-End",
        "smarts": "[CH3][NH][CH]([CH2][CH]1[CH2][NH][C](=[NH])[NH]1)[C](=[O])[OH]",
    },
    {
        "name": "C16:0-OH(3)",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH]([OH])[CH2][C](=[O])[OH]",
    },
    {
        "name": "NFo-Cl-Trp",
        "smarts": "[O]=[CH][NH][CH]([CH2][C]1=[N][c]2[cH][c]([Cl])[cH][cH][c]2[CH2]1)[C](=[O])[OH]",
    },
    {
        "name": "C10:0-Me(2.4)-oxo(9)",
        "smarts": "[CH3][C](=[O])[CH2][CH2][CH2][CH2][CH]([CH3])[CH2][CH]([CH3])[C](=[O])[OH]",
    },
    {"name": "Hpa", "smarts": "[O]=[C]([OH])[c]1[n][cH][cH][cH][c]1[OH]"},
    {
        "name": "iC11:0-OH(3)",
        "smarts": "[CH3][CH]([CH3])[CH2][CH2][CH2][CH2][CH2][CH]([OH])[CH2][C](=[O])[OH]",
    },
    {"name": "Aib", "smarts": "[CH3][C]([CH3])([NH2])[C](=[O])[OH]"},
    {
        "name": "C8:0:1(7)-Me(2)-OH(3)",
        "smarts": "[CH]#[C][CH2][CH2][CH2][CH]([OH])[CH]([CH3])[C](=[O])[OH]",
    },
    {"name": "NFo-Dab", "smarts": "[NH2][CH2][CH2][CH]([NH][CH]=[O])[C](=[O])[OH]"},
    {"name": "4oxo-Hpr", "smarts": "[O]=[C]1[CH2][CH2][NH][CH]([C](=[O])[OH])[CH2]1"},
    {
        "name": "Phe-Thz",
        "smarts": "[NH2][CH]([CH2][c]1[cH][cH][cH][cH][cH]1)[c]1[n][c]([C](=[O])[OH])[cH][s]1",
    },
    {
        "name": "NAc-1Me-Trp",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][c]1[cH][n]([CH3])[c]2[cH][cH][cH][cH][c]12)[C](=[O])[OH]",
    },
    {
        "name": "NMe-bMe-Br-Phe",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[CH]([CH3])[c]1[cH][cH][c]([Br])[cH][cH]1",
    },
    {
        "name": "NMe-OMe-TyrC",
        "smarts": "[CH3][NH][CH]([CH2][c]1[cH][cH][c]([O][CH3])[cH][cH]1)[C]([NH2])=[O]",
    },
    {
        "name": "NFo-ChrI",
        "smarts": "[NH2][C]1=[CH][c]2[cH][c]([OH])[c]([OH])[cH][c]2[N]2[CH2][CH2][CH]([C](=[O])[OH])[N]([CH]=[O])[CH]12",
    },
    {
        "name": "b5OH-Cap",
        "smarts": "[NH]=[C]1[NH][CH]([OH])[CH2][C]([NH2])([CH2][C](=[O])[OH])[NH]1",
    },
    {
        "name": "NFo-bOMe-Asp",
        "smarts": "[CH3][O][C](=[O])[CH2][CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "NAc-4oxo-5Me-Pro",
        "smarts": "[CH3][C](=[O])[N]1[CH]([CH3])[C](=[O])[CH2][CH]1[C](=[O])[OH]",
    },
    {"name": "Glu", "smarts": "[NH2][CH]([CH2][CH2][C](=[O])[OH])[C](=[O])[OH]"},
    {
        "name": "bOH-Br-Phe",
        "smarts": "[NH2][CH]([C](=[O])[OH])[CH]([OH])[c]1[cH][cH][c]([Br])[cH][cH]1",
    },
    {
        "name": "3Me-Glu",
        "smarts": "[CH3][CH]([CH2][C](=[O])[OH])[CH]([NH2])[C](=[O])[OH]",
    },
    {
        "name": "MeO-Glu",
        "smarts": "[CH3][O][C](=[O])[CH2][CH2][CH]([NH2])[C](=[O])[OH]",
    },
    {
        "name": "e-Tyr",
        "smarts": "[NH2][CH]([CH2][CH2][C](=[O])[OH])[CH2][c]1[cH][cH][c]([OH])[cH][cH]1",
    },
    {
        "name": "C16:0-NH2(3)",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH]([NH2])[CH2][C](=[O])[OH]",
    },
    {
        "name": "NFo-Cl-CONH2-Trp",
        "smarts": "[NH2][C](=[O])[n]1[cH][c]([CH2][CH]([NH][CH]=[O])[C](=[O])[OH])[c]2[cH][cH][c]([Cl])[cH][c]21",
    },
    {"name": "bOH-Asp", "smarts": "[NH2][C]([OH])([CH2][C](=[O])[OH])[C](=[O])[OH]"},
    {"name": "Eta", "smarts": "[NH2][CH2][CH2][OH]"},
    {
        "name": "NMe-bMe-Asn",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[CH]([CH3])[C]([NH2])=[O]",
    },
    {"name": "PT", "smarts": "[CH3][P](=[O])([OH])[CH2][CH2][CH]([NH2])[C](=[O])[OH]"},
    {
        "name": "NFo-Aca",
        "smarts": "[O]=[CH][NH][CH]([CH2][CH]1[CH2][CH2][C](=[O])[CH]2[O][CH]12)[C](=[O])[OH]",
    },
    {
        "name": "bOH-NMe-Phe",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[CH]([OH])[c]1[cH][cH][cH][cH][cH]1",
    },
    {"name": "bDab", "smarts": "[NH2][CH2][CH]([NH2])[CH2][C](=[O])[OH]"},
    {
        "name": "NFo-Orn",
        "smarts": "[NH2][CH2][CH2][CH2][CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "NAc-Asp",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][C](=[O])[OH])[C](=[O])[OH]",
    },
    {
        "name": "NAc-OMe-Thr",
        "smarts": "[CH3][O][CH]([CH3])[CH]([NH][C]([CH3])=[O])[C](=[O])[OH]",
    },
    {
        "name": "Cl2-NMe-dhLeu",
        "smarts": "[CH3][NH][C](=[CH][CH]([CH3])[CH]([Cl])[Cl])[C](=[O])[OH]",
    },
    {"name": "Hpg", "smarts": "[NH2][CH]([C](=[O])[OH])[c]1[cH][cH][c]([OH])[cH][cH]1"},
    {
        "name": "NMe-bMe-Asp",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[CH]([CH3])[C](=[O])[OH]",
    },
    {"name": "bGlu", "smarts": "[NH2][CH]([CH2][C](=[O])[OH])[CH2][C](=[O])[OH]"},
    {"name": "bU-dAla", "smarts": "[NH2][C](=[O])[NH][CH]=[C]([NH2])[C](=[O])[OH]"},
    {
        "name": "pTrp",
        "smarts": "[O]=[C]([OH])[CH]1[CH2][C]2([OH])[c]3[cH][cH][cH][cH][c]3[NH][CH]2[NH]1",
    },
    {
        "name": "NFo-Phe",
        "smarts": "[O]=[CH][NH][CH]([CH2][c]1[cH][cH][cH][cH][cH]1)[C](=[O])[OH]",
    },
    {
        "name": "b3Me-Phe",
        "smarts": "[CH3][c]1[cH][cH][cH][c]([CH]([NH2])[CH2][C](=[O])[OH])[cH]1",
    },
    {
        "name": "bC10:0-NH2(2)-Ep(9)-oxo(8)",
        "smarts": "[NH2][CH]([CH2][CH2][CH2][CH2][C](=[O])[CH]1[CH2][O]1)[CH2][C](=[O])[OH]",
    },
    {
        "name": "NMe-Fo-OH-Orn",
        "smarts": "[CH3][NH][CH]([CH2][CH2][CH2][N]([OH])[CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "bbOH-Cl-Tyr",
        "smarts": "[NH2][C]([OH])([CH2][C](=[O])[OH])[c]1[cH][cH][c]([OH])[c]([Cl])[cH]1",
    },
    {
        "name": "NFo-1Me-Trp",
        "smarts": "[CH3][n]1[cH][c]([CH2][CH]([NH][CH]=[O])[C](=[O])[OH])[c]2[cH][cH][cH][cH][c]21",
    },
    {
        "name": "bbOH-Br-Phe",
        "smarts": "[NH2][C]([OH])([CH2][C](=[O])[OH])[c]1[cH][cH][c]([Br])[cH][cH]1",
    },
    {
        "name": "NFo-bMe-Ile",
        "smarts": "[CH3][CH2][C]([CH3])([CH3])[CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "bbNMe-NMe-Trp",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[C]([CH3])([CH3])[c]1[cH][n]([CH3])[c]2[cH][cH][cH][cH][c]12",
    },
    {
        "name": "C10:0-OH(9)-NH2(2)-oxo(8)",
        "smarts": "[CH3][CH]([OH])[C](=[O])[CH2][CH2][CH2][CH2][CH2][CH]([NH2])[C](=[O])[OH]",
    },
    {
        "name": "C18:1(9)",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH]=[CH][CH2][CH2][CH2][CH2][CH2][CH2][CH2][C](=[O])[OH]",
    },
    {
        "name": "NFo-Glu",
        "smarts": "[O]=[CH][NH][CH]([CH2][CH2][C](=[O])[OH])[C](=[O])[OH]",
    },
    {
        "name": "iC15:0-OH(3)",
        "smarts": "[CH3][CH]([CH3])[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH]([OH])[CH2][C](=[O])[OH]",
    },
    {"name": "Cys", "smarts": "[NH2][CH]([CH2][SH])[C](=[O])[OH]"},
    {
        "name": "NMe-N1-COOH-bhTrp",
        "smarts": "[CH3][NH][CH]([CH2][CH2][CH2][c]1[cH][n]([C](=[O])[OH])[c]2[cH][cH][cH][cH][c]12)[C](=[O])[OH]",
    },
    {
        "name": "NAc-PT",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][CH2][P]([CH3])(=[O])[OH])[C](=[O])[OH]",
    },
    {
        "name": "gOH-NMe-Val",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[CH]([CH3])[CH2][OH]",
    },
    {
        "name": "NMe-OMe-Thr",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[CH]([CH3])[O][CH3]",
    },
    {
        "name": "bTyr",
        "smarts": "[NH2][CH]([CH2][C](=[O])[OH])[c]1[cH][cH][c]([OH])[cH][cH]1",
    },
    {"name": "OH-Orn", "smarts": "[NH2][CH]([CH2][CH2][CH2][NH][OH])[C](=[O])[OH]"},
    {"name": "NFo-Asn", "smarts": "[NH2][C](=[O])[CH2][CH]([NH][CH]=[O])[C](=[O])[OH]"},
    {
        "name": "C4:0-OH(2.3)-Cl(4)",
        "smarts": "[O]=[C]([OH])[CH]([OH])[CH]([OH])[CH2][Cl]",
    },
    {
        "name": "NFo-OMe-Thr",
        "smarts": "[CH3][O][CH]([CH3])[CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {"name": "k-Leu", "smarts": "[CH3][CH]([CH3])[CH2][C](=[O])[C](=[O])[OH]"},
    {"name": "bMet", "smarts": "[CH3][S][CH2][CH]([NH2])[CH2][C](=[O])[OH]"},
    {"name": "aThr/Thr", "smarts": "[CH3][CH]([OH])[CH]([NH2])[C](=[O])[OH]"},
    {"name": "Rha", "smarts": "[CH3][CH]1[O][CH]([OH])[CH]([OH])[CH]([OH])[CH]1[OH]"},
    {"name": "Ala", "smarts": "[CH3][CH]([NH2])[C](=[O])[OH]"},
    {
        "name": "CFA",
        "smarts": "[CH3][CH2][CH]1[CH]=[C]([C](=[O])[OH])[CH]2[CH2][CH2][C](=[O])[CH]2[CH2]1",
    },
    {
        "name": "NFo-4Me-Pro",
        "smarts": "[CH3][CH]1[CH2][CH]([C](=[O])[OH])[N]([CH]=[O])[CH2]1",
    },
    {
        "name": "NAc-4OH-Thr",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[CH]([OH])[CH2][OH]",
    },
    {
        "name": "NMe-pTrp",
        "smarts": "[CH3][N]1[CH]([C](=[O])[OH])[CH2][C]2([OH])[c]3[cH][cH][cH][cH][c]3[NH][CH]12",
    },
    {"name": "Lys", "smarts": "[NH2][CH2][CH2][CH2][CH2][CH]([NH2])[C](=[O])[OH]"},
    {
        "name": "NAc-bOH-Val",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[C]([CH3])([CH3])[OH]",
    },
    {
        "name": "NMe-Arg",
        "smarts": "[CH3][NH][CH]([CH2][CH2][CH2][N]=[C]([NH2])[NH2])[C](=[O])[OH]",
    },
    {"name": "ProC", "smarts": "[NH2][C](=[O])[CH]1[CH2][CH2][CH2][NH]1"},
    {
        "name": "Apv",
        "smarts": "[NH2][CH]([CH2][CH2][CH2][c]1[cH][cH][cH][cH][cH]1)[C](=[O])[OH]",
    },
    {
        "name": "NFo-Apv",
        "smarts": "[O]=[CH][NH][CH]([CH2][CH2][CH2][c]1[cH][cH][cH][cH][cH]1)[C](=[O])[OH]",
    },
    {
        "name": "NAc-Ph-Ser",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[CH]([OH])[c]1[cH][cH][cH][cH][cH]1",
    },
    {
        "name": "NFo-3Me-Phe",
        "smarts": "[CH3][c]1[cH][cH][cH][c]([CH2][CH]([NH][CH]=[O])[C](=[O])[OH])[cH]1",
    },
    {
        "name": "NMe-C10:0-NH2(2)-Ep(9)-oxo(8)",
        "smarts": "[CH3][NH][CH]([CH2][CH2][CH2][CH2][CH2][C](=[O])[CH]1[CH2][O]1)[C](=[O])[OH]",
    },
    {
        "name": "NAc-N1-COOH-bhTrp",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][CH2][CH2][c]1[cH][n]([C](=[O])[OH])[c]2[cH][cH][cH][cH][c]12)[C](=[O])[OH]",
    },
    {"name": "bHis", "smarts": "[NH2][CH]([CH2][C](=[O])[OH])[c]1[cH][n][cH][nH]1"},
    {
        "name": "C9:0-Me(2)-OH(3)",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH]([OH])[CH]([CH3])[C](=[O])[OH]",
    },
    {
        "name": "OMe-Trp",
        "smarts": "[CH3][O][n]1[cH][c]([CH2][CH]([NH2])[C](=[O])[OH])[c]2[cH][cH][cH][cH][c]21",
    },
    {
        "name": "NMe-4Cl-Thr",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[CH]([OH])[CH2][Cl]",
    },
    {
        "name": "NFo-Bmt",
        "smarts": "[CH3][CH]=[CH][CH2][CH]([CH3])[CH]([OH])[CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "C12:3(7.9.11)-Me(6)-OH(2.4.5)-NH2(3)-Ph(12)",
        "smarts": "[CH3][CH]([CH]=[CH][CH]=[CH][CH]=[CH][c]1[cH][cH][cH][cH][cH]1)[CH]([OH])[CH]([OH])[CH]([NH2])[CH]([OH])[C](=[O])[OH]",
    },
    {"name": "NMe-Asp", "smarts": "[CH3][NH][CH]([CH2][C](=[O])[OH])[C](=[O])[OH]"},
    {"name": "bOH-Asn", "smarts": "[NH2][C](=[O])[C]([NH2])([OH])[CH2][C](=[O])[OH]"},
    {
        "name": "bC10:0-OH(8)-NH2(2)",
        "smarts": "[CH3][CH2][CH]([OH])[CH2][CH2][CH2][CH2][CH]([NH2])[CH2][C](=[O])[OH]",
    },
    {"name": "bbMe-Asn", "smarts": "[CH3][C]([NH2])([CH2][C](=[O])[OH])[C]([NH2])=[O]"},
    {
        "name": "Argal",
        "smarts": "[NH2][C]([NH2])=[N][CH2][CH2][CH2][CH]([NH2])[CH]=[O]",
    },
    {
        "name": "Gal/Man/bGal/Glc",
        "smarts": "[OH][CH2][CH]1[O][CH]([OH])[CH]([OH])[CH]([OH])[CH]1[OH]",
    },
    {
        "name": "bBmt",
        "smarts": "[CH3][CH]=[CH][CH2][CH]([CH3])[C]([NH2])([OH])[CH2][C](=[O])[OH]",
    },
    {
        "name": "NFo-Cl2-Pro",
        "smarts": "[O]=[CH][N]1[CH2][CH]([Cl])[CH]([Cl])[CH]1[C](=[O])[OH]",
    },
    {"name": "4oxo-Pro", "smarts": "[O]=[C]1[CH2][NH][CH]([C](=[O])[OH])[CH2]1"},
    {
        "name": "aC13:1(3)",
        "smarts": "[CH3][CH2][CH]([CH3])[CH2][CH2][CH2][CH2][CH2][CH]=[CH][CH2][C](=[O])[OH]",
    },
    {"name": "Bz", "smarts": "[O]=[C]([OH])[c]1[cH][cH][cH][cH][cH]1"},
    {"name": "bOrn", "smarts": "[NH2][CH2][CH2][CH]([NH2])[CH2][C](=[O])[OH]"},
    {
        "name": "NAc-Thr",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[CH]([CH3])[OH]",
    },
    {
        "name": "NMe-Me-AOA",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH]([NH][CH3])[C](=[O])[OH]",
    },
    {
        "name": "NFo-3Me-Glu",
        "smarts": "[CH3][CH]([CH2][C](=[O])[OH])[CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "NMe-CysA",
        "smarts": "[CH3][NH][CH]([CH2][S](=[O])(=[O])[OH])[C](=[O])[OH]",
    },
    {
        "name": "NMe-PO-Asn",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[CH]([O][P](=[O])([OH])[OH])[C]([NH2])=[O]",
    },
    {"name": "O-Met", "smarts": "[CH3][S](=[O])[CH2][CH2][CH]([NH2])[C](=[O])[OH]"},
    {
        "name": "NMe-4OH-Pro",
        "smarts": "[CH3][N]1[CH2][CH]([OH])[CH2][CH]1[C](=[O])[OH]",
    },
    {
        "name": "NFo-Cit",
        "smarts": "[NH2][C](=[O])[NH][CH2][CH2][CH2][CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "Ist",
        "smarts": "[CH3][CH2][CH]([CH3])[CH]([NH2])[CH]([OH])[CH2][C](=[O])[OH]",
    },
    {
        "name": "C10:2(7.9)-OH(2.4.5)-NH2(3)-ePh(10)",
        "smarts": "[CH3][CH2][O][c]1[cH][cH][c]([CH]=[CH][CH]=[CH][CH2][CH]([OH])[CH]([OH])[CH]([NH2])[CH]([OH])[C](=[O])[OH])[cH][cH]1",
    },
    {
        "name": "NAc-Br-OH-Trp",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][c]1[c]([Br])[nH][c]2[cH][cH][c]([OH])[cH][c]12)[C](=[O])[OH]",
    },
    {
        "name": "Cl-Tyr",
        "smarts": "[NH2][CH]([CH2][c]1[cH][cH][c]([OH])[c]([Cl])[cH]1)[C](=[O])[OH]",
    },
    {
        "name": "NFo-bOH-Tyr",
        "smarts": "[O]=[CH][NH][CH]([C](=[O])[OH])[CH]([OH])[c]1[cH][cH][c]([OH])[cH][cH]1",
    },
    {
        "name": "NMe-Hpg",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[c]1[cH][cH][c]([OH])[cH][cH]1",
    },
    {
        "name": "NAc-O2-Met",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][CH2][S]([CH3])(=[O])=[O])[C](=[O])[OH]",
    },
    {"name": "bPT", "smarts": "[CH3][P](=[O])([OH])[CH2][CH]([NH2])[CH2][C](=[O])[OH]"},
    {"name": "NMe-Cys", "smarts": "[CH3][NH][CH]([CH2][SH])[C](=[O])[OH]"},
    {"name": "Ac-Ser", "smarts": "[CH3][C](=[O])[O][CH2][CH]([NH2])[C](=[O])[OH]"},
    {
        "name": "Cl-Hpg",
        "smarts": "[NH2][CH]([C](=[O])[OH])[c]1[cH][cH][c]([OH])[c]([Cl])[cH]1",
    },
    {
        "name": "NAc-Br-Phe",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][c]1[cH][cH][c]([Br])[cH][cH]1)[C](=[O])[OH]",
    },
    {"name": "NMe-Met", "smarts": "[CH3][NH][CH]([CH2][CH2][S][CH3])[C](=[O])[OH]"},
    {
        "name": "NAc-Ahad",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][CH]([OH])[CH2][C](=[O])[OH])[C](=[O])[OH]",
    },
    {
        "name": "NAc-3Me-Phe",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][c]1[cH][cH][cH][c]([CH3])[cH]1)[C](=[O])[OH]",
    },
    {
        "name": "Agdha",
        "smarts": "[NH]=[C]([NH2])[NH][CH2][CH2][CH2][CH]([NH2])[CH]([OH])[CH]([OH])[C](=[O])[OH]",
    },
    {
        "name": "NMe-COOH-Trp",
        "smarts": "[CH3][NH][CH]([CH2][c]1[c]([C](=[O])[OH])[nH][c]2[cH][cH][cH][cH][c]12)[C](=[O])[OH]",
    },
    {
        "name": "NFo-Ph-Ser",
        "smarts": "[O]=[CH][NH][CH]([C](=[O])[OH])[CH]([OH])[c]1[cH][cH][cH][cH][cH]1",
    },
    {
        "name": "Hty",
        "smarts": "[NH2][CH]([CH2][CH2][c]1[cH][cH][c]([OH])[cH][cH]1)[C](=[O])[OH]",
    },
    {
        "name": "C8:2(5.7)-Me(6)-OH(4)-NH2(3)-Ph(8)",
        "smarts": "[CH3][C]([CH]=[CH][c]1[cH][cH][cH][cH][cH]1)=[CH][CH]([OH])[CH]([NH2])[CH2][C](=[O])[OH]",
    },
    {
        "name": "Kyn",
        "smarts": "[NH2][c]1[cH][cH][cH][cH][c]1[C](=[O])[CH2][CH]([NH2])[C](=[O])[OH]",
    },
    {"name": "aMe-Cys", "smarts": "[CH3][C]([NH2])([CH2][SH])[C](=[O])[OH]"},
    {
        "name": "Nst",
        "smarts": "[CH3][CH]([CH3])[CH2][CH]([NH2])[CH]([OH])[C](=[O])[OH]",
    },
    {"name": "Azd", "smarts": "[O]=[C]([OH])[CH]1[NH][CH]1[C](=[O])[OH]"},
    {
        "name": "NMe-Hil",
        "smarts": "[CH3][CH2][CH]([CH3])[CH2][CH]([NH][CH3])[C](=[O])[OH]",
    },
    {"name": "NFo-Dpr", "smarts": "[NH2][CH2][CH]([NH][CH]=[O])[C](=[O])[OH]"},
    {
        "name": "NFo-Br-Trp",
        "smarts": "[O]=[CH][NH][CH]([CH2][c]1[cH][nH][c]2[cH][cH][c]([Br])[cH][c]12)[C](=[O])[OH]",
    },
    {
        "name": "bMe-AOA",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH]([NH2])[CH2][C](=[O])[OH]",
    },
    {
        "name": "NAc-Aca",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][CH]1[CH2][CH2][C](=[O])[CH]2[O][CH]12)[C](=[O])[OH]",
    },
    {"name": "bMe-Asn", "smarts": "[CH3][CH]([C]([NH2])=[O])[CH]([NH2])[C](=[O])[OH]"},
    {
        "name": "C12:1(5)-OH(3)",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH]=[CH][CH2][CH]([OH])[CH2][C](=[O])[OH]",
    },
    {
        "name": "4OH-Ph-Lac",
        "smarts": "[O]=[C]([OH])[CH]([OH])[CH2][c]1[cH][cH][c]([OH])[cH][cH]1",
    },
    {
        "name": "NMe-3Me-Pro",
        "smarts": "[CH3][CH]1[CH2][CH2][N]([CH3])[CH]1[C](=[O])[OH]",
    },
    {
        "name": "NAc-3Me-Hty",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[CH]([CH3])[CH2][c]1[cH][cH][c]([OH])[cH][cH]1",
    },
    {
        "name": "NFo-bOH-Br-Phe",
        "smarts": "[O]=[CH][NH][CH]([C](=[O])[OH])[CH]([OH])[c]1[cH][cH][c]([Br])[cH][cH]1",
    },
    {
        "name": "NAc-Daz",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][CH2][CH2][CH]([NH2])[CH]([OH])[CH2][C](=[O])[OH])[C](=[O])[OH]",
    },
    {
        "name": "NAc-Tyr",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][c]1[cH][cH][c]([OH])[cH][cH]1)[C](=[O])[OH]",
    },
    {"name": "NAc-Ser", "smarts": "[CH3][C](=[O])[NH][CH]([CH2][OH])[C](=[O])[OH]"},
    {
        "name": "b3Me-Hty",
        "smarts": "[CH3][C]([NH2])([CH2][C](=[O])[OH])[CH2][c]1[cH][cH][c]([OH])[cH][cH]1",
    },
    {"name": "Phe", "smarts": "[NH2][CH]([CH2][c]1[cH][cH][cH][cH][cH]1)[C](=[O])[OH]"},
    {
        "name": "C7:0-OH(3)",
        "smarts": "[CH3][CH2][CH2][CH2][CH]([OH])[CH2][C](=[O])[OH]",
    },
    {"name": "Val", "smarts": "[CH3][CH]([CH3])[CH]([NH2])[C](=[O])[OH]"},
    {"name": "dPyr", "smarts": "[NH2][CH]1[CH2][C](=[O])[NH][C]1=[CH][C](=[O])[OH]"},
    {
        "name": "iC17:0-OH(3)",
        "smarts": "[CH3][CH]([CH3])[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH]([OH])[CH2][C](=[O])[OH]",
    },
    {
        "name": "bOMe-Tyr",
        "smarts": "[CH3][O][CH]([c]1[cH][cH][c]([OH])[cH][cH]1)[CH]([NH2])[C](=[O])[OH]",
    },
    {
        "name": "NFo-4oxo-Hpr",
        "smarts": "[O]=[CH][N]1[CH2][CH2][C](=[O])[CH2][CH]1[C](=[O])[OH]",
    },
    {
        "name": "NMe-3Me-Hty",
        "smarts": "[CH3][NH][CH]([C](=[O])[OH])[CH]([CH3])[CH2][c]1[cH][cH][c]([OH])[cH][cH]1",
    },
    {
        "name": "bPh-Ser",
        "smarts": "[NH2][C]([OH])([CH2][C](=[O])[OH])[c]1[cH][cH][cH][cH][cH]1",
    },
    {
        "name": "OSu-Hmp",
        "smarts": "[CH3][CH2][CH]([CH3])[CH]([O][S](=[O])(=[O])[OH])[C](=[O])[OH]",
    },
    {"name": "NSpd", "smarts": "[NH2][CH2][CH2][CH2][NH][CH2][CH2][CH2][NH2]"},
    {
        "name": "NMe-Trp",
        "smarts": "[CH3][NH][CH]([CH2][c]1[cH][nH][c]2[cH][cH][cH][cH][c]12)[C](=[O])[OH]",
    },
    {
        "name": "4oxo-5Me-Pro",
        "smarts": "[CH3][CH]1[NH][CH]([C](=[O])[OH])[CH2][C]1=[O]",
    },
    {
        "name": "NAc-Glu",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][CH2][C](=[O])[OH])[C](=[O])[OH]",
    },
    {
        "name": "NAc-bOMe-Asp",
        "smarts": "[CH3][O][C](=[O])[CH2][CH]([NH][C]([CH3])=[O])[C](=[O])[OH]",
    },
    {
        "name": "ChrP",
        "smarts": "[NH2][C]1=[CH][c]2[cH][c]([OH])[c]([OH])[cH][c]2[N]2[CH]([C](=[O])[OH])[CH2][CH2][NH][CH]12",
    },
    {
        "name": "C6:0-OH(3.5)-NH2(4)",
        "smarts": "[CH3][CH]([OH])[CH]([NH2])[CH]([OH])[CH2][C](=[O])[OH]",
    },
    {"name": "Cl2-Pro", "smarts": "[O]=[C]([OH])[CH]1[NH][CH2][CH]([Cl])[CH]1[Cl]"},
    {
        "name": "aC17:0-OH(3)",
        "smarts": "[CH3][CH2][CH]([CH3])[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH]([OH])[CH2][C](=[O])[OH]",
    },
    {
        "name": "NFo-O-Met",
        "smarts": "[CH3][S](=[O])[CH2][CH2][CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "NFo-bOMe-Tyr",
        "smarts": "[CH3][O][CH]([c]1[cH][cH][c]([OH])[cH][cH]1)[CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {"name": "C6:0", "smarts": "[CH3][CH2][CH2][CH2][CH2][C](=[O])[OH]"},
    {
        "name": "NAc-Cl-CONH2-Trp",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][c]1[cH][n]([C]([NH2])=[O])[c]2[cH][c]([Cl])[cH][cH][c]12)[C](=[O])[OH]",
    },
    {"name": "C8:0:1(7)", "smarts": "[CH]#[C][CH2][CH2][CH2][CH2][CH2][C](=[O])[OH]"},
    {
        "name": "NFo-Gln",
        "smarts": "[NH2][C](=[O])[CH2][CH2][CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "NAc-Gln",
        "smarts": "[CH3][C](=[O])[NH][CH]([CH2][CH2][C]([NH2])=[O])[C](=[O])[OH]",
    },
    {
        "name": "NFo-His",
        "smarts": "[O]=[CH][NH][CH]([CH2][c]1[cH][n][cH][nH]1)[C](=[O])[OH]",
    },
    {
        "name": "NFo-F-ph-Gly",
        "smarts": "[O]=[CH][NH][CH]([C](=[O])[OH])[c]1[cH][cH][c]([F])[cH][cH]1",
    },
    {
        "name": "NFo-PT",
        "smarts": "[CH3][P](=[O])([OH])[CH2][CH2][CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "NAc-bOH-Cl-Tyr",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[CH]([OH])[c]1[cH][cH][c]([OH])[c]([Cl])[cH]1",
    },
    {
        "name": "NAc-OH-His",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[CH]([OH])[c]1[cH][n][cH][nH]1",
    },
    {
        "name": "C10:0-NH2(2)-Ep(9)-oxo(8)",
        "smarts": "[NH2][CH]([CH2][CH2][CH2][CH2][CH2][C](=[O])[CH]1[CH2][O]1)[C](=[O])[OH]",
    },
    {
        "name": "Ahad",
        "smarts": "[NH2][CH]([CH2][CH]([OH])[CH2][C](=[O])[OH])[C](=[O])[OH]",
    },
    {
        "name": "C14:0-NH2(3)",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH]([NH2])[CH2][C](=[O])[OH]",
    },
    {
        "name": "NAc-3OH-5Me-Pro",
        "smarts": "[CH3][C](=[O])[N]1[CH]([CH3])[CH2][CH]([OH])[CH]1[C](=[O])[OH]",
    },
    {
        "name": "iC14:0",
        "smarts": "[CH3][CH]([CH3])[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][C](=[O])[OH]",
    },
    {
        "name": "C10:0-Me(4)-OH(3)",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH]([CH3])[CH]([OH])[CH2][C](=[O])[OH]",
    },
    {"name": "NMe-Dpr", "smarts": "[CH3][NH][CH]([CH2][NH2])[C](=[O])[OH]"},
    {
        "name": "NFo-bMe-Gln",
        "smarts": "[CH3][CH]([CH2][C]([NH2])=[O])[CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {"name": "3OH-Pro", "smarts": "[O]=[C]([OH])[CH]1[NH][CH2][CH2][CH]1[OH]"},
    {
        "name": "U4oxo-Van",
        "smarts": "[CH3][CH]1[O][CH]([OH])[CH2][C]2([CH3])[NH][C](=[O])[NH][C]12[OH]",
    },
    {
        "name": "Ibu",
        "smarts": "[CH3][CH]([NH2])[C](=[O])[C]([CH3])([CH3])[C](=[O])[OH]",
    },
    {
        "name": "NFo-MeO-Glu",
        "smarts": "[CH3][O][C](=[O])[CH2][CH2][CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "NAc-bMe-Phe",
        "smarts": "[CH3][C](=[O])[NH][CH]([C](=[O])[OH])[CH]([CH3])[c]1[cH][cH][cH][cH][cH]1",
    },
    {
        "name": "DHMDA",
        "smarts": "[CH3][CH2][CH2][CH]([CH3])[CH]([OH])[CH]([CH3])[CH]([CH]=[CH][CH]=[CH][CH2][C](=[O])[OH])[O][CH3]",
    },
    {
        "name": "NFo-Lys",
        "smarts": "[NH2][CH2][CH2][CH2][CH2][CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {"name": "NFo-Ala", "smarts": "[CH3][CH]([NH][CH]=[O])[C](=[O])[OH]"},
    {"name": "iC5:0-OH(2.3)", "smarts": "[CH3][C]([CH3])([OH])[CH]([OH])[C](=[O])[OH]"},
    {"name": "NMe-Ser", "smarts": "[CH3][NH][CH]([CH2][OH])[C](=[O])[OH]"},
    {"name": "NMe-Pro", "smarts": "[CH3][N]1[CH2][CH2][CH2][CH]1[C](=[O])[OH]"},
    {
        "name": "bHar",
        "smarts": "[NH2][C]([NH2])=[N][CH2][CH2][CH2][CH]([NH2])[CH2][C](=[O])[OH]",
    },
    {"name": "NMe-Hpr", "smarts": "[CH3][N]1[CH2][CH2][CH2][CH2][CH]1[C](=[O])[OH]"},
    {
        "name": "Daz",
        "smarts": "[NH2][CH]([CH2][CH2][CH2][CH]([NH2])[CH]([OH])[CH2][C](=[O])[OH])[C](=[O])[OH]",
    },
    {
        "name": "bCl-Trp",
        "smarts": "[NH2][CH]([CH2][C](=[O])[OH])[C]1=[N][c]2[cH][c]([Cl])[cH][cH][c]2[CH2]1",
    },
    {
        "name": "NAc-3Me-Pro",
        "smarts": "[CH3][C](=[O])[N]1[CH2][CH2][CH]([CH3])[CH]1[C](=[O])[OH]",
    },
    {
        "name": "NFo-Et-Nva",
        "smarts": "[CH3][CH2][CH]([CH2][CH3])[CH]([NH][CH]=[O])[C](=[O])[OH]",
    },
    {
        "name": "C13:2(t4.t6)-OH(2.3)",
        "smarts": "[CH3][CH2][CH2][CH2][CH2][CH2][CH]=[CH][CH]=[CH][CH]([OH])[CH]([OH])[C](=[O])[OH]",
    },
    {
        "name": "NMe-Gln",
        "smarts": "[CH3][NH][CH]([CH2][CH2][C]([NH2])=[O])[C](=[O])[OH]",
    },
]
