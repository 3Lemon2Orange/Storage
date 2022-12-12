// Print Month.
#include <iostream>

using namespace std;

// determine whether a year is a leap year
bool is_leap_year(int year)
{
    return (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
}

// return number of days in a month
int get_num_month_days(int year, int month)
{
    if (month == 4 || month == 6 || month == 9 || month == 11) return 30;
    else if (month == 2) {
        if (is_leap_year(year)) return 29;
        else return 28;
    }
    else return 31;
}

// return the anchor day of a century
int get_anchor_century(int year)
{
    int century;
    century = year / 100;
    switch (century % 4) {
        case 0:
            return 2;
        case 1:
            return 0;
        case 2:
            return 5;
        case 3:
            return 3;
    }
}

// return the weekday after offset
int offset_weekday(int weekday, int offset)
{
    return ((weekday + offset) % 7 + 7) % 7;
}

// return the anchor day of a year (the weekday of doomsdays)
int get_anchor_year(int year, int anchor_century)
{
    int yy, a, b, c;
    yy = year % 100;
    a = yy / 12;
    b = yy % 12;
    c = b / 4;
    return offset_weekday(anchor_century, a + b + c);
}

// return the doomsday in a month
int get_doomsday(int year, int month)
{
    if (month == 1) {
        if (is_leap_year(year)) return 4;
        else return 3;
    }
    else if (month == 2) {
        if (is_leap_year(year)) return 1;
        else return 7;
    }
    else if (month == 3 || month == 5 || month == 7)
        return month + 4;
    else if (month == 9 || month == 11)
        return month - 4;
    else return month;
}

// return the weekday of a date
int get_weekday(int year, int month, int day)
{
    int anchor_century, anchor_year, doomsday;
    anchor_century = get_anchor_century(year);
    anchor_year = get_anchor_year(year, anchor_century);
    doomsday = get_doomsday(year, month);
    return offset_weekday(anchor_year, day - doomsday);
}

// print the calendar of the month containing the date
void print_calendar(int year, int month, int day)
{
    // print title
    cout << "==========================================================" << endl;
    cout << "Doomsday Calendar, in memorial of John Conway (1937-2020)." << endl;
    cout << "==========================================================" << endl;
    cout << "\t\t\t" <<year << "-" << month << endl;
    cout << "Sun\tMon\tTus\tWed\tThu\tFri\tSat" << endl;
    // print blank spaces
    int weekday = get_weekday(year, month, 1);
    for (int i = 0; i < weekday; i++) {
        cout << "\t";
    }
    // print dates
    int num_month_days = get_num_month_days(year, month);
    for (int date = 1; date <= num_month_days; date++) {
        // print date
        cout << date;
        if (date == day) cout << "*";
        // formating
        if (date < num_month_days) {
            if (weekday == 6) cout << endl;
            else cout << "\t";
        }
        else cout << endl;
        // update weekday
        weekday = offset_weekday(weekday, 1);
    }
    cout << "==========================================================" << endl;
}

int main()
{
    int year, month, day;
    cout << "Please input the date ('year month day'): ";
    cin >> year >> month >> day;
    print_calendar(year, month, day);
    return 0;
}
