-- Schema for temperaturecycling

create table temperature (
    id          integer primary key autoincrement not null,
    added_time  datetime default current_timestamp,
    base_temp   int
);