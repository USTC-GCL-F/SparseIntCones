#include "CommonFunction.h"

CommonFunction::CommonFunction()
{

}

vector<vector<int>> CommonFunction::enumeration(int var_num,int diff)
{
	vector<vector<int>> guass_bonet;
	switch (var_num)
	{
		case 1:guass_bonet = enum_var1(var_num, diff); break;
		case 2:guass_bonet = enum_var2(var_num, diff); break;
		case 3:guass_bonet = enum_var3(var_num, diff); break;
		case 4:guass_bonet = enum_var4(var_num, diff); break;
		case 5:guass_bonet = enum_var5(var_num, diff); break;
		case 6:guass_bonet = enum_var6(var_num, diff); break;
		case 7:guass_bonet = enum_var7(var_num, diff); break;
		case 8:guass_bonet = enum_var8(var_num, diff); break;
		case 9:guass_bonet = enum_var9(var_num, diff); break;
		case 10:guass_bonet = enum_var10(var_num, diff); break;
		default:cout << endl;
				cout << "变量个数超过10个" << endl;
				abort();
				exit(EXIT_FAILURE);;
	}
	return guass_bonet;
}

vector<vector<int>> CommonFunction::enumeration_high_genus(int var_num, int diff)
{
	vector<vector<int>> guass_bonet;
	switch (var_num)
	{
	case 1:guass_bonet = enum_var1(var_num, diff); break;
	case 2:guass_bonet = enum_var2(var_num, diff); break;
	case 3:guass_bonet = enum_var3(var_num, diff); break;
	case 4:guass_bonet = enum_var4(var_num, diff); break;
	case 5:guass_bonet = enum_var5(var_num, diff); break;
	case 6:guass_bonet = enum_var6(var_num, diff); break;
	case 7:guass_bonet = enum_var7(var_num, diff); break;
	case 8:guass_bonet = enum_var8(var_num, diff); break;
	case 9:guass_bonet = enum_var9(var_num, diff); break;
	case 10:guass_bonet = enum_var10(var_num, diff); break;
	case 11:guass_bonet = enum_var11(var_num, diff); break;
	case 12:guass_bonet = enum_var12(var_num, diff); break;
	case 13:guass_bonet = enum_var13(var_num, diff); break;
	case 14:guass_bonet = enum_var14(var_num, diff); break;
	default:cout << endl;
		cout << "变量个数超过14个" << endl;
		abort();
		exit(EXIT_FAILURE);;
	}
	return guass_bonet;
}



vector<vector<int>> CommonFunction::enum_var1(int var_num, int diff)
{
	vector<int> num = { -1,0,1 };
	vector<vector<int>> guass_bonet;
	vector<int> sub_vec;
	for (int i = 0; i < 3; i++)
	{
		sub_vec = { num[i] };
		int sum_vec = accumulate(sub_vec.begin(), sub_vec.end(), 0);
		if (sum_vec== diff)
		{			
			guass_bonet.push_back(sub_vec);			
		}
		sub_vec.clear();	
	}
	return guass_bonet;
}

vector<vector<int>> CommonFunction::enum_var2(int var_num, int diff)
{
	vector<int> num = { -1,0,1 };
	vector<vector<int>> guass_bonet;
	vector<int> sub_vec;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			sub_vec = { num[i] ,num[j] };
			int sum_vec = accumulate(sub_vec.begin(), sub_vec.end(), 0);
			if (sum_vec == diff)
			{
				guass_bonet.push_back(sub_vec);
			}
			sub_vec.clear();
		}
	}
	return guass_bonet;
}

vector<vector<int>> CommonFunction::enum_var3(int var_num, int diff)
{
	vector<int> num = { -1,0,1 };
	vector<vector<int>> guass_bonet;
	vector<int> sub_vec;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				sub_vec = { num[i] ,num[j] ,num[k]};
				int sum_vec = accumulate(sub_vec.begin(), sub_vec.end(), 0);
				if (sum_vec == diff)
				{
					guass_bonet.push_back(sub_vec);
				}
				sub_vec.clear();
			}
		}
	}
	return guass_bonet;	
}

vector<vector<int>> CommonFunction::enum_var4(int var_num, int diff)
{
	vector<int> num = { -1,0,1 };
	vector<vector<int>> guass_bonet;
	vector<int> sub_vec;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				for (int l = 0; l < 3; l++)
				{
					sub_vec = { num[i] ,num[j] ,num[k],num[l]};
					int sum_vec = accumulate(sub_vec.begin(), sub_vec.end(), 0);
					if (sum_vec == diff)
					{
						guass_bonet.push_back(sub_vec);
					}
					sub_vec.clear();
				}
			}
		}
	}
	return guass_bonet;
}

vector<vector<int>> CommonFunction::enum_var5(int var_num, int diff)
{
	vector<int> num = { -1,0,1 };
	vector<vector<int>> guass_bonet;
	vector<int> sub_vec;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				for (int l = 0; l < 3; l++)
				{
					for (int z = 0; z < 3; z++)
					{
						sub_vec = { num[i] ,num[j] ,num[k],num[l],num[z]};
						int sum_vec = accumulate(sub_vec.begin(), sub_vec.end(), 0);
						if (sum_vec == diff)
						{
							guass_bonet.push_back(sub_vec);
						}
						sub_vec.clear();
					}
				}
			}
		}
	}
	return guass_bonet;
}

vector<vector<int>> CommonFunction::enum_var6(int var_num, int diff)
{
	vector<int> num = { -1,0,1 };
	vector<vector<int>> guass_bonet;
	vector<int> sub_vec;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				for (int l = 0; l < 3; l++)
				{
					for (int z = 0; z < 3; z++)
					{
						for (int q = 0; q < 3; q++)
						{
							sub_vec = { num[i] ,num[j] ,num[k],num[l],num[z] ,num[q]};
							int sum_vec = accumulate(sub_vec.begin(), sub_vec.end(), 0);
							if (sum_vec == diff)
							{
								guass_bonet.push_back(sub_vec);
							}
							sub_vec.clear();
						}
					}
				}
			}
		}
	}
	return guass_bonet;
}

vector<vector<int>> CommonFunction::enum_var7(int var_num, int diff)
{
	vector<int> num = { -1,0,1 };
	vector<vector<int>> guass_bonet;
	vector<int> sub_vec;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				for (int l = 0; l < 3; l++)
				{
					for (int z = 0; z < 3; z++)
					{
						for (int q = 0; q < 3; q++)
						{
							for (int q7 = 0; q7 < 3; q7++)
							{
								sub_vec = { num[i] ,num[j] ,num[k],num[l],num[z] ,num[q] ,num[q7]};
								int sum_vec = accumulate(sub_vec.begin(), sub_vec.end(), 0);
								if (sum_vec == diff)
								{
									guass_bonet.push_back(sub_vec);
								}
								sub_vec.clear();
							}
						}
					}
				}
			}
		}
	}
	return guass_bonet;
}

vector<vector<int>> CommonFunction::enum_var8(int var_num, int diff)
{
	vector<int> num = { -1,0,1 };
	vector<vector<int>> guass_bonet;
	vector<int> sub_vec;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				for (int l = 0; l < 3; l++)
				{
					for (int z = 0; z < 3; z++)
					{
						for (int q = 0; q < 3; q++)
						{
							for (int q7 = 0; q7 < 3; q7++)
							{
								for (int q8 = 0; q8 < 3; q8++)
								{
									sub_vec = { num[i] ,num[j] ,num[k],num[l],num[z] ,num[q] ,num[q7] ,num[q8]};
									int sum_vec = accumulate(sub_vec.begin(), sub_vec.end(), 0);
									if (sum_vec == diff)
									{
										guass_bonet.push_back(sub_vec);
									}
									sub_vec.clear();
								}
							}
						}
					}
				}
			}
		}
	}
	return guass_bonet;
}

vector<vector<int>> CommonFunction::enum_var9(int var_num, int diff)
{
	vector<int> num = { -1,0,1 };
	vector<vector<int>> guass_bonet;
	vector<int> sub_vec;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				for (int l = 0; l < 3; l++)
				{
					for (int z = 0; z < 3; z++)
					{
						for (int q = 0; q < 3; q++)
						{
							for (int q7 = 0; q7 < 3; q7++)
							{
								for (int q8 = 0; q8 < 3; q8++)
								{
									for (int q9 = 0; q9 < 3; q9++)
									{
										sub_vec = { num[i] ,num[j] ,num[k],num[l],num[z] ,num[q] ,num[q7] ,num[q8],num[q9]};
										int sum_vec = accumulate(sub_vec.begin(), sub_vec.end(), 0);
										if (sum_vec == diff)
										{
											guass_bonet.push_back(sub_vec);
										}
										sub_vec.clear();
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return guass_bonet;
}

vector<vector<int>> CommonFunction::enum_var10(int var_num, int diff)
{
	vector<int> num = { -1,0,1 };
	vector<vector<int>> guass_bonet;
	vector<int> sub_vec;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				for (int l = 0; l < 3; l++)
				{
					for (int z = 0; z < 3; z++)
					{
						for (int q = 0; q < 3; q++)
						{
							for (int q7 = 0; q7 < 3; q7++)
							{
								for (int q8 = 0; q8 < 3; q8++)
								{
									for (int q9 = 0; q9 < 3; q9++)
									{
										for (int q10 = 0; q10 < 3; q10++)
										{
											sub_vec = { num[i] ,num[j] ,num[k],num[l],num[z] ,num[q] ,num[q7] ,num[q8],num[q9] ,num[q10]};
											int sum_vec = accumulate(sub_vec.begin(), sub_vec.end(), 0);
											if (sum_vec == diff)
											{
												guass_bonet.push_back(sub_vec);
											}
											sub_vec.clear();
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return guass_bonet;
}

vector<vector<int>> CommonFunction::enum_var11(int var_num, int diff)
{
	vector<int> num = { -1,0,1 };
	vector<vector<int>> guass_bonet;
	vector<int> sub_vec;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				for (int l = 0; l < 3; l++)
				{
					for (int z = 0; z < 3; z++)
					{
						for (int q = 0; q < 3; q++)
						{
							for (int q7 = 0; q7 < 3; q7++)
							{
								for (int q8 = 0; q8 < 3; q8++)
								{
									for (int q9 = 0; q9 < 3; q9++)
									{
										for (int q10 = 0; q10 < 3; q10++)
										{
											for (int q11 = 0; q11 < 3; q11++)
											{
												sub_vec = { num[i] ,num[j] ,num[k],num[l],num[z] ,num[q] ,num[q7] ,num[q8],num[q9] ,num[q10] ,num[q11] };
												int sum_vec = accumulate(sub_vec.begin(), sub_vec.end(), 0);
												if (sum_vec == diff)
												{
													guass_bonet.push_back(sub_vec);
												}
												sub_vec.clear();
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return guass_bonet;
}

vector<vector<int>> CommonFunction::enum_var12(int var_num, int diff)
{
	vector<int> num = { -1,0,1 };
	vector<vector<int>> guass_bonet;
	vector<int> sub_vec;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				for (int l = 0; l < 3; l++)
				{
					for (int z = 0; z < 3; z++)
					{
						for (int q = 0; q < 3; q++)
						{
							for (int q7 = 0; q7 < 3; q7++)
							{
								for (int q8 = 0; q8 < 3; q8++)
								{
									for (int q9 = 0; q9 < 3; q9++)
									{
										for (int q10 = 0; q10 < 3; q10++)
										{
											for (int q11 = 0; q11 < 3; q11++)
											{
												for (int q12 = 0; q12 < 3; q12++)
												{	
													sub_vec = { num[i] ,num[j] ,num[k],num[l],num[z] ,num[q] ,num[q7] ,num[q8],num[q9] ,num[q10] ,num[q11],num[q12] };
													int sum_vec = accumulate(sub_vec.begin(), sub_vec.end(), 0);
													if (sum_vec == diff)
													{
														guass_bonet.push_back(sub_vec);
													}
													sub_vec.clear();
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return guass_bonet;
}

vector<vector<int>> CommonFunction::enum_var13(int var_num, int diff)
{
	vector<int> num = { -1,0,1 };
	vector<vector<int>> guass_bonet;
	vector<int> sub_vec;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				for (int l = 0; l < 3; l++)
				{
					for (int z = 0; z < 3; z++)
					{
						for (int q = 0; q < 3; q++)
						{
							for (int q7 = 0; q7 < 3; q7++)
							{
								for (int q8 = 0; q8 < 3; q8++)
								{
									for (int q9 = 0; q9 < 3; q9++)
									{
										for (int q10 = 0; q10 < 3; q10++)
										{
											for (int q11 = 0; q11 < 3; q11++)
											{
												for (int q12 = 0; q12 < 3; q12++)
												{
													for (int q13 = 0; q13 < 3; q13++)
													{
														sub_vec = { num[i] ,num[j] ,num[k],num[l],num[z] ,num[q] ,num[q7] ,num[q8],num[q9] ,num[q10] ,num[q11],num[q12] ,num[q13] };
														int sum_vec = accumulate(sub_vec.begin(), sub_vec.end(), 0);
														if (sum_vec == diff)
														{
															guass_bonet.push_back(sub_vec);
														}
														sub_vec.clear();
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return guass_bonet;
}

vector<vector<int>> CommonFunction::enum_var14(int var_num, int diff)
{
	vector<int> num = { -1,0,1 };
	vector<vector<int>> guass_bonet;
	vector<int> sub_vec;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				for (int l = 0; l < 3; l++)
				{
					for (int z = 0; z < 3; z++)
					{
						for (int q = 0; q < 3; q++)
						{
							for (int q7 = 0; q7 < 3; q7++)
							{
								for (int q8 = 0; q8 < 3; q8++)
								{
									for (int q9 = 0; q9 < 3; q9++)
									{
										for (int q10 = 0; q10 < 3; q10++)
										{
											for (int q11 = 0; q11 < 3; q11++)
											{
												for (int q12 = 0; q12 < 3; q12++)
												{
													for (int q13 = 0; q13 < 3; q13++)
													{
														for (int q14 = 0; q14 < 3; q14++)
														{
															sub_vec = { num[i] ,num[j] ,num[k],num[l],num[z] ,num[q] ,num[q7] ,num[q8],num[q9] ,num[q10] ,num[q11],num[q12] ,num[q13] ,num[q14] };
															int sum_vec = accumulate(sub_vec.begin(), sub_vec.end(), 0);
															if (sum_vec == diff)
															{
																guass_bonet.push_back(sub_vec);
															}
															sub_vec.clear();
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return guass_bonet;
}























